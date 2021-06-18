#ifndef MultiSUPERPIXEL_HIERARCHY_MEX_HPP
#define MultiSUPERPIXEL_HIERARCHY_MEX_HPP

#include <algorithm>
#include <iostream>
#include <string.h>
//#include <cmath>
#include <math.h>       /* sqrt */
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

template<class T>
void radixSortLSD(int *I, int *Itemp, T *J, T *Jtemp, int length) // radixSortLSD(m_dist,m_dtemp,m_minarc,m_arctmp,m_vexnum);
{
	int *beginI  = I,     *endI  =     I+length;
	int *beginIt = Itemp, *endIt = Itemp+length;
	T   *beginJ  = J,     *endJ  =     J+length;
	T   *beginJt = Jtemp, *endJt = Jtemp+length;
	for (int shift = 0; shift < 32; shift += 8) 
	{
		size_t count[0x100] = {};
		for (int *p = beginI; p != endI; p++)
			count[(*p >> shift) & 0xFF]++;
		int *bucketI[0x100], *qI = beginIt;
		T   *bucketJ[0x100], *qJ = beginJt;
		for (int i = 0; i < 0x100; ++i)
		{
			bucketI[i] = qI; bucketJ[i] = qJ;
			qI += count[i]; qJ += count[i];
		}
		int *p = beginI;
		T   *q = beginJ;
		for (; p != endI; p++, q++)
		{
			int idx = (*p >> shift) & 0xFF;
			*bucketI[idx]++ = *p;
			*bucketJ[idx]++ = *q;
		}
		std::swap(beginI, beginIt);
		std::swap(beginJ, beginJt);
		std::swap(endI, endIt);
		std::swap(endJ, endJt);
	}
}

struct ArcBox
{
	int u, v;
	ArcBox *ulink, *vlink;
	unsigned short confidence;
	//unsigned char   ismin;
	//int length;
	//ArcBox() : u(0), v(0), ulink(0), vlink(0), confidence(0) {} //
};

class MultiSuperpixelHierarchy
{

private:
	ArcBox *m_arcmem;
	ArcBox **m_arcptr;
	ArcBox **m_frtarc;
	ArcBox **m_minarc;
	ArcBox **m_arctmp;

	int *m_vex, *m_parent, *m_label, *m_size, *m_newvex;
	int *m_dist, m_vexmax, m_arcmax, m_vexnum, m_arcnum;
	int m_h, m_w, m_connect, m_iter, m_iterSwitch;
	int m_maxDistColor, m_chnl, m_selPerc, m_maxSize;
	int m_joinedpixels;
	double m_mean;
	bool m_euclidean, m_meanColor, joinedContinue;

	unsigned short *m_color, *m_edge;
	int *m_treeu, *m_treev, *m_dtemp;
	int m_treeSize, m_regionnum; 

	int *m_temp;

public:
	MultiSuperpixelHierarchy()
	{
		m_arcmem = NULL;
		m_arcptr = NULL;
		m_frtarc = NULL;
		m_minarc = NULL;
		m_arctmp = NULL;
		m_color  = NULL;
		m_temp   = NULL;
		m_edge   = NULL;

	}

	~MultiSuperpixelHierarchy() { clean(); }

	void init(int h, int w, int chnl, int connect = 24, int iterSwitch = 0, int SelPerc = 100, int maxDistColor = 5, int maxSize = 10000, bool euclidean = true, bool meanColor = true)
	{
		//	SH.init(h,w,chnl,connect,iterSwitch,SelPerc,SparseJoin,maxSize,euclidean);

		m_h = h; m_w = w; 
		m_connect = connect; m_iterSwitch = iterSwitch; m_chnl = chnl; m_selPerc = SelPerc; m_maxDistColor = maxDistColor; m_maxSize = maxSize; m_euclidean=euclidean, m_meanColor=meanColor;
		m_vexmax = m_vexnum = h*w;
		m_arcmax = computeEdge(h,w,connect); // max number of arcs
		m_arcnum = 0;
		m_arcmem = new ArcBox [m_arcmax];
		m_arcptr = new ArcBox*[m_arcmax];
		m_minarc = new ArcBox*[m_vexnum]; 
		m_arctmp = new ArcBox*[m_vexnum];
		m_frtarc = new ArcBox*[m_vexnum]; 
		memset(m_frtarc,0,sizeof(ArcBox*)*m_vexnum);
		
		m_color  = new unsigned short[m_vexnum*m_chnl];
		m_edge = new unsigned short[m_vexnum];
		m_temp   = new int[8*m_vexnum]; 
		m_vex    = m_temp+0*m_vexnum; 
		m_newvex = m_vex;
		m_size   = m_temp+1*m_vexnum; 
		m_parent = m_temp+2*m_vexnum; 
		m_label  = m_temp+3*m_vexnum;
		m_dtemp  = m_temp+4*m_vexnum;  
		m_dist   = m_temp+5*m_vexnum; 
		m_treeu  = m_temp+6*m_vexnum;  
		m_treev  = m_temp+7*m_vexnum; 

		for (int i=0; i<m_vexnum; ++i) { m_vex[i] = m_parent[i] = m_label[i] = i; m_size[i] = 1;}
		for (int i=0; i<m_arcmax; ++i) { m_arcptr[i] = &m_arcmem[i]; } 
		m_regionnum = m_vexnum; m_treeSize = 0; m_joinedpixels = 0; joinedContinue  = true;

	}

	void buildTree(unsigned short *img, unsigned short *edge = NULL)
	{

		createVexNchnl(img);
		createVexedge(edge);
		buildGraph(edge);

		m_iter = 0;
		//int maxDistColor = 0;
        // mexPrintf("Lets find nearestNegihbor");
		// mexEvalString("drawnow;");
		
		while (m_vexnum>5) // threshold is not met
		{
			++m_iter;
			nearestNeighbor();
//			if (m_SparseJoin)
//				growRegionJoin();
//			else
			growRegion();			
			merge();

			// mexPrintf("buildtree() - m_iter", m_iter);
			// mexEvalString("drawnow;");
		}
	}

	bool insertArc(int u, int v, unsigned short confidence = 0)
	{
		//if (isNeighbor(i, j)) return false;
		//if (i > j) std::swap(i, j); // let i be the small vertex

		ArcBox *arc  = m_arcptr[m_arcnum++];
		arc->u = u;
		arc->v = v;
		arc->ulink = m_frtarc[u];
		arc->vlink = m_frtarc[v];
		arc->confidence = confidence;
		m_frtarc[u] = m_frtarc[v] = arc;

		// mexPrintf("arc-u = %i", arc->u);
		// mexEvalString("drawnow;");

		return true;
	}

	ArcBox *findArc(ArcBox *arc, int u, int v)
	{
		while (arc)
		{
			if (arc->v == v)
				return arc;
			else if (arc->v == u)
				arc = arc->vlink;
			else
				arc = arc->ulink;
		}
		return NULL;
	}

	void nearestNeighbor()
	{
		// mexPrintf("Find nearest \n");
		// mexEvalString("drawnow;");
		double sum = 0;
		for (int i=0; i<m_vexnum; ++i)
		{
			sum += m_size[i];
		}
		m_mean=sum/m_vexnum;	

		for (int i=0; i<m_vexnum; ++i) m_dist[i] = 2147483645;
		for (int n=0; n<m_arcnum; ++n)
		{
			ArcBox *arc  = m_arcptr[n];
			int u = arc->u;
			int v = arc->v;
			int lu = m_label[u];
			int lv = m_label[v];
			
			unsigned short *uc = m_color+m_chnl*u;
			unsigned short *vc = m_color+m_chnl*v;
			
			double distColor;
			// mexPrintf("uc vc = %i.%i", *uc , *vc);
			// mexEvalString("drawnow;");
            
			//if (m_euclidean)
			//{
				//	[1] color term - euclidean dist. 
				double dc[m_chnl];
				for (int ii=0; ii<m_chnl; ++ii){
					dc[ii] = uc[ii]-vc[ii];	
				}
				distColor = 0;
				for (int ii=0; ii<m_chnl; ++ii){
					distColor += dc[ii]*dc[ii];	
				}

				

				//distColor = distColor;
				// error is from 0 to 1.000.000. m_dist of 50.000 is 5% of error.
				//distColor = sqrt(distColor);
				// mexPrintf("Euclideandist-%f, ",distColor);
				// mexEvalString("drawnow;");
				// distColor /= m_chnl;
				//				
				// distColor /= sqrt(m_chnl)*65535;
				distColor = sqrt(distColor)*(1000000/(sqrt(m_chnl)*65535));///65536/m_chnl; 
		
				//if (m_iter > m_iterSwitch) distColor /= arc->confidence;
				// mexPrintf("dist-%f pu-%i pv-%i,",distColor, u, v);
				// mexEvalString("drawnow;");

				distColor *= ((m_size[u]+m_size[v])/(m_mean*2));//*((m_size[u]+m_size[v])/(m_mean*2));						

				
				// mexPrintf("%g,", distColor);
				// mexEvalString("drawnow;");
				
			//}
			//else
			//{
				//	[1] color term - cosine dist.			
			//	distColor = cosineDist(uc,vc);	
			//	distColor = (distColor*1000000);///65536/m_chnl; 				
			//	if (m_iter > m_iterSwitch) distColor /= arc->confidence;						

			//}			
							

			// [2] boundary term
			// int distEdge = arc->confidence;		
			// Regularize distances using  auxiliary images.
			// distColor /= m_edge[i];						
			
			
			// Regularize distances to avoid little superpixels.
			//if (m_iter > m_iterSwitch) distColor *= exp((double)(((m_size[u]+m_size[v])-m_maxSize)/m_maxSize));
			
			//mexPrintf("%e,", exp((double)(((m_size[u]+m_size[v])-m_maxSize)/m_maxSize)));
			//mexEvalString("drawnow;");
 				

			int dist = (int)distColor;


			// // Regularize distances to avoid very big superpixels.
			// if ((m_size[u]>m_maxSize)||(m_size[v]>m_maxSize) )
			// {
		 	// 	dist = 100000000; 
			// }

			if (dist < m_dist[lu]) { m_dist[lu] = dist; m_minarc[lu] = arc; }
			if (dist < m_dist[lv]) { m_dist[lv] = dist; m_minarc[lv] = arc; }

			//mexPrintf("%i ", m_dist[lu]);
			//mexEvalString("drawnow;");
		}

	}

	void growRegion()
	{

		// radix sort. Hace sorting de los arcs segun su dist!
		radixSortLSD(m_dist,m_dtemp,m_minarc,m_arctmp,m_vexnum);
		mexPrintf("Number of superpixels = %i \n",m_vexnum);
		mexEvalString("drawnow;");
		if (m_vexnum>100 & m_iter>1 & joinedContinue)
		{	

			// mexPrintf("std+mean = %i.",(int)std+mean);
			// mexEvalString("drawnow;");
			for (int i=0; i<ceil((m_vexnum*m_selPerc)/100); ++i)
			{				
				ArcBox *arc = m_minarc[i];
				
				int pu = findset(arc->u);
				int pv = findset(arc->v);				
				// int lu = m_label[findset(arc->u)];
				// int lv = m_label[findset(arc->v)];
				//int lu = findset(arc->u);
				//int lv = findset(arc->v);
				// mexPrintf("dist-%i pu-%i pv-%i,",m_dist[i], pu, pv);
				// mexPrintf("distNormalized-%f, lu-%i, lv-%i",(m_dist[i]/((m_size[lu]+m_size[lv])/(m_mean*2))),lu,lv);
				// mexEvalString("drawnow;");
				if (pu == pv) continue;				

				//if (((m_dist[i]/((m_size[lu]+m_size[lv])/(m_mean*2)))>(m_maxDistColor*10000)))
				// if ((m_dist[i]*(m_size[pu]+m_size[pv])>(m_maxDistColor*10000)))
				// {					
					// if (i>5){
						// m_joinedpixels = m_treeSize;
						// mexPrintf("Sppxl Threshold: %i, mdist=%i,", m_joinedpixels,m_dist[i]);
				   	 	// mexEvalString("drawnow;");
						// break; 
					// }
						
					// if (i<5){
						// joinedContinue=false;
						// m_joinedpixels = m_treeSize;
						// mexPrintf("Sppxl Threshold: %i, mdist=%i,", m_joinedpixels,m_dist[i]);
				   	 	// mexEvalString("drawnow;");
						// break; 
					// }


						
				// }
				
				// if (m_dist[i]>(int)(mean) & m_joinedpixels==0)
				// {	
				// 	// mexPrintf("mdist,  i: %i,%i\n", m_dist[i],i);
				//     // mexEvalString("drawnow;");
				// 	break;
				// }

				m_treeu[m_treeSize]   = pu;
				m_treev[m_treeSize++] = pv;
				// m_countVex[pu]++;
				// m_countVex[pv]++;
				// mexPrintf("treeSize=%i, ",m_treeSize);
				// mexEvalString("drawnow;");
				if (pu < pv) m_parent[pv] = pu;
				else         m_parent[pu] = pv;
			}
			
			// update parent and label			
			// // Do randomSample.	
			// int p[m_vexnum-1] = {};
			// for (int i=0; i<m_vexnum; ++i) {
			// 	(p[i]=i);
			// 	// mexPrintf("%i ", p[i]);
			// 	// mexEvalString("drawnow;");
			// }
			// //mexPrintf("%i ", p[m_vexnum]);
			// for (int i=m_vexnum-1; i>0; --i)
			// {
			// 	//get swap index
			// 	int j = rand()%i;
			// 	//swap p[i] with p[j]
			// 	int temp = p[i];
			// 	p[i] = p[j];
			// 	p[j] = temp;				
			// }

			// mexPrintf("RandomSample done\n");
		    // mexEvalString("drawnow;");
			m_regionnum = 0;			
			for (int i=0; i<m_vexnum; ++i)
			{
				// mexPrintf("%i ", p[i]);
		    	// mexEvalString("drawnow;");	
				int u = m_vex[i];
				int pu = findset(u);
				
				// mexPrintf("u = %i, pu = %i.\n", u, pu);
				// mexEvalString("drawnow;");				

				if (pu == u)
				{
					m_label[u] = m_regionnum;
					// create a new component
					m_vex[m_regionnum] = u;
					m_frtarc[m_regionnum] = NULL; //Eliminate arc
					++m_regionnum;
				}
				else
				{
					m_label[u] = m_label[pu];
					
					
					//if (m_meanColor)
					//{
						// update color
						int s1 = m_size[u]; int s2 = m_size[pu]; int s = s1+s2;
						m_size[pu] = s;
						// Mean color of both.
						for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = (m_color[m_chnl*u+ii]*s1 + m_color[m_chnl*pu+ii]*s2)/s;
						//mexPrintf("u-%i, m_color(1)=%i", u, m_color[m_chnl*u+0]);
						//mexEvalString("drawnow;");
					//}
					//else
					//{
						// // Max density of both
						// if (m_edge[u]>m_edge[pu])
						// {
						// 	for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = m_color[m_chnl*u+ii];
						// 	m_edge[pu]=m_edge[u];
						// }
						
						// Weighted color by density and size of pixels.
					//	int s1 = m_size[u]; int s2 = m_size[pu]; int s = s1+s2;
					//	m_size[pu] = s;
					//	int d1 = m_edge[u]; int d2 = m_edge[pu]; int d = d1+d2; 
						//m_edge[pu] = (s1*e1+s2*e2)/s; 
					//	m_edge[pu] = (((m_edge[u]*s1+m_edge[pu]*s2)/(s1+s2))+((m_edge[u]*d1+m_edge[pu]*d2)/(d1+d2)))/2;
						
						// for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = (m_color[m_chnl*u+ii]*s1 + m_color[m_chnl*pu+ii]*s2)/s;
						// for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = (m_color[m_chnl*u+ii]*e1 + m_color[m_chnl*pu+ii]*e2)/e;

					//	for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = (((m_color[m_chnl*u+ii]*s1+m_color[m_chnl*pu+ii]*s2)/(s1+s2))+((m_color[m_chnl*u+ii]*d1+m_color[m_chnl*pu+ii]*d2)/(d1+d2)))/2;

						// mexPrintf("e1-%i, m_edge = %i, m_color(1)=%i\n", u, m_edge[u], m_color[m_chnl*u+0]);
						// mexEvalString("drawnow;");

						//int maxU = m_color[m_chnl*u];
						//int maxPU = m_color[m_chnl*pu];
						// Calculating the min of each superpixel
						//for (int ii=1; ii<m_chnl; ++ii) 
						//{
						//	if (m_color[m_chnl*u+ii] > maxU) maxU = m_color[m_chnl*u+ii];
						//	if (m_color[m_chnl*pu+ii] > maxPU) maxPU = m_color[m_chnl*pu+ii];
						//}					
						//if (maxU>maxPU) for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = m_color[m_chnl*u+ii];

					//}
					
				}
				
			}
			//m_vex = m_newvex;
			m_vexnum = m_regionnum;
			
			mexPrintf("vex_num %i", m_vexnum);
			mexEvalString("drawnow;");
		}
		else
		{
			for (int i=0; i<m_vexnum; ++i)
			{
				// mexPrintf("dist = %i, i = %i.\n", m_dist[i],i);
				// mexEvalString("drawnow;");
				
				ArcBox *arc = m_minarc[i];
				// mexPrintf("Arc %i\n",arc->u);
				// mexEvalString("drawnow;");
				
				int pu = findset(arc->u);
				int pv = findset(arc->v);
				// mexPrintf("dist = %i, pu = %i, arcpu = %i.\n", m_dist[i], pu, arc->u);
				// mexEvalString("drawnow;");
				if (pu == pv) continue;

				m_treeu[m_treeSize]   = pu;
				m_treev[m_treeSize++] = pv;
				// m_countVex[pu]++;
				// m_countVex[pv]++;
				// mexPrintf("pu=%i, counting=%i \n",pu ,m_countVex[pu]);
				// mexEvalString("drawnow;");
				if (pu < pv) m_parent[pv] = pu;
				else         m_parent[pu] = pv;
			}
			
		
			// update parent and label
			m_regionnum = 0;
			for (int i=0; i<m_vexnum; ++i)
			{
				int u = m_vex[i];
				int pu = findset(u);
				
				// mexPrintf("u = %i, pu = %i.\n", u, pu);
				// mexEvalString("drawnow;");

				if (pu == u)
				{
					m_label[u] = m_regionnum;
					// create a new component
					m_vex[m_regionnum] = u;
					m_frtarc[m_regionnum] = NULL; //Eliminate arc
					++m_regionnum;
				}
				else
				{
					m_label[u] = m_label[pu];
					// update color
					int s1 = m_size[u]; int s2 = m_size[pu]; int s = s1+s2;
					m_size[pu] = s;
					// Mean color of both.
					for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = (m_color[m_chnl*u+ii]*s1 + m_color[m_chnl*pu+ii]*s2)/s;
				}
				
			}

			m_vexnum = m_regionnum;

		}
	}

	void merge()
	{
		// update connection
		// mexPrintf("Update Connection\n");
		// mexEvalString("drawnow;");
		int arcold = m_arcnum;
		m_arcnum = 0;
		for (int n=0; n<arcold; ++n)
		{
			// mexPrintf("%i,", m_dist[n]);
			// mexEvalString("drawnow;");

			// mexPrintf("arcs %i of %i.\n", n,arcold);
			// mexEvalString("drawnow;");

			ArcBox *arc = m_arcptr[n];
			int  u = arc->u;
			int  v = arc->v;
			int pu = m_parent[u];
			int pv = m_parent[v];
			// mexPrintf("u %i, v %i, pu %i, pv %i.\n", u, v, pu,pv);
			// mexEvalString("drawnow;");
			if (pu == pv)
				continue;

			if (pu > pv) { std::swap(u,v); std::swap(pu,pv); }
			int lu = m_label[u]; int lv = m_label[v];
			ArcBox *&arcu = m_frtarc[lu];
			ArcBox *&arcv = m_frtarc[lv];
			ArcBox * find = findArc(arcu, pu, pv);
			if (find)
			{
				// find->confidence = (m_size[pu]+m_size[u]);
				if (find->confidence > arc->confidence)
					// find->confidence = arc->confidence;
					find->confidence = m_edge[pu];
			}
			else
			{
				//assert(arc->ismin == 0);
				arc->u = pu;
				arc->v = pv;
				arc->ulink = arcu;
				arc->vlink = arcv;
				//arc->confidence = (m_size[pu]+m_size[u]);; 
				m_arcptr[m_arcnum++] = arcu = arcv = arc;
			}
			// mexPrintf("arcs=%i-%i. confidence=%i\n", arc->u, arc->v, arc->confidence);
			// mexEvalString("drawnow;");
		}
		// mexPrintf("Finished Update connection\n");
		// mexEvalString("drawnow;");
	}

	int findset(int i)
	{
		// mexPrintf("Parent[i], %i\n", m_parent[i]);
		// mexEvalString("drawnow;");
		int p = m_parent[i];
		if (i != p)
		{
			m_parent[i] = findset(p);
		}
		return m_parent[i];
	}

	int *getLabel()  { return m_label;     }
	int *getParent() { return m_parent;    }
	int *getTreeU()  { return m_treeu;     }
	int *getTreeV()  { return m_treev;     }
	int  getRegion() { return m_regionnum; }
	int getJoinedPixels() { return m_joinedpixels; }

	int *getLabel(int N)
	{
		if (N < 1 || N > m_vexmax)
		{
			printf("error");
			exit(1);
		}

		int end   = m_vexmax-N;
		int begin = m_vexmax-m_regionnum;
		if (m_regionnum < N)
		{
			for (int i=0; i<m_vexmax; ++i) m_parent[i] = i;
			begin = 0;
		}

		for (int i=begin; i<end; ++i)
		{
			int u  = m_treeu[i];
			int v  = m_treev[i];
			int pu = findset(u);
			int pv = findset(v);
			if (pu < pv)
				m_parent[pv] = pu;
			else
				m_parent[pu] = pv;
		}

		m_regionnum = 0;
		for (int i=0; i<m_vexmax; ++i)
		{
			int p = findset(i);
			if (i == p)
				m_label[i] = m_regionnum++;
			else
				m_label[i] = m_label[p];
		}

		return m_label;
	}

private:

	void clean()
	{
		delete [] m_arcptr;
		delete [] m_arcmem;
		delete [] m_frtarc;
		delete [] m_minarc;
		delete [] m_arctmp;
		delete [] m_color;
		delete [] m_edge;
		delete [] m_temp;
	}

	int computeEdge(int h, int w, int connect)
	{
		if (h<=2 || w<=2)
			return 0;
		else if (connect == 4)
			return (h-1)*w + (w-1)*h;
		else if (connect == 8)
			return (h-1)*w + (w-1)*h + (h-1)*(w-1)*2;	
		else if (connect == 24)
			return (h-1)*w + (w-1)*h + (h-1)*(w-1)*2 + (h-2)*(w)*3 + (h)*(w-2)*3 + (h-2)*(w-2)*2;	
		else if (connect == 48)
		 	return (h-1)*w + (w-1)*h + (h-1)*(w-1)*2;	
	}


	void createVexNchnl(unsigned short *img)
	{
		unsigned short  *ptrr  = img;
		unsigned short  *ptrg  = ptrr+m_vexnum;
		unsigned short  *ptre  = ptrg;
		unsigned short *ptrd  = m_color;
		while (ptrr != ptre)
		{
			for (int n=1; n<m_chnl; ++n)
			{
				*(ptrd+n) = *(ptrr+m_vexnum*n);					
			}
			*ptrd = *ptrr++;
			ptrd += m_chnl;
		}
	}

	void createVexedge(unsigned short *edge)
	{
		unsigned short  *ptrr  = edge;
		unsigned short  *ptrg  = ptrr+m_vexnum;
		unsigned short  *ptre  = ptrg;
		unsigned short  *ptrd  = m_edge;
		while (ptrr != ptre)
		{
			*ptrd = *ptrr++;
			ptrd += 1;
		}
	}

	void buildGraph(unsigned short *edge)
	{		
		// x direction (calculo de edges que van hacia x) hacia abajo.
		for (int y=0; y<m_h; ++y)
		{
			int id0 = y;
			int idt = id0+m_h;
			for (int x=0; x<m_w-1; ++x)
			{
				unsigned short boundary = min(edge[id0], edge[idt]); //el minimo de los extremos es el boundary.
				// mexPrintf("xpabajo = %i.\n",boundary );
				// mexEvalString("drawnow;");
				//boundary = boundary + 1;
				//boundary = max(boundary,1); // poner a uno si es menor.
				insertArc(id0,idt,boundary);
				id0 = idt; idt += m_h;
			}
		}

		// y direction (hacia la derecha)
		for (int x=0; x<m_w; ++x)
		{
			int id0 = x*m_h;
			int idt = id0+1;
			for (int y=0; y<m_h-1; ++y)
			{
				unsigned short boundary = min(edge[id0], edge[idt]);
				//boundary = boundary + 1;
				boundary = max(boundary,1);
				insertArc(id0,idt,boundary);
				id0 = idt++;
			}
		}

		if (m_connect >= 8)
		{
			for (int y0=0; y0<m_h-1; ++y0)
			{
				int yt = y0+1;
				for (int x0=0; x0<m_w-1; ++x0)
				{
					int xt = x0+1;
					int id0 = x0*m_h+y0;
					int idt = xt*m_h+yt;
					unsigned short boundary = min(edge[id0], edge[idt]);
					//boundary = boundary + 1;
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					id0 = xt*m_h+y0;
					idt = x0*m_h+yt;
					boundary = min(edge[id0], edge[idt]);
					//boundary = boundary + 1;
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
				}
			}
		}

		if (m_connect >= 24)
		{	
			//Siendo x=id0  x-1 
			//            8/|\2
			//            76543
			for (int y0=0; y0<m_h-2; ++y0)
			{
				int yt1 = y0+2; // 1,2,3 
				int yt2 = y0+1; // 4,8
				int yt3 = y0; // 5,6,7
				for (int x0=0; x0<m_w-2; ++x0)
				{
					int xt1 = x0; // 1
					int xt2 = x0+1; // 2,6
					int xt3 = x0+2; // 3,4,5,7,8
					int id0 = x0*m_h+y0;
					// 1.
					int idt = xt1*m_h+yt1;
					unsigned short boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 2.
					idt = xt2*m_h+yt1;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 3.
					idt = xt3*m_h+yt1;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 4.
					idt = xt3*m_h+yt2;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 5.
					idt = xt3*m_h+yt3;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 6.
					id0 = xt2*m_h+y0;
					idt = x0*m_h+yt3;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 7.
					id0 = xt3*m_h+y0;
					idt = x0*m_h+yt3;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
					// 8.
					id0 = xt3*m_h+y0;
					idt = x0*m_h+yt2;
					boundary = min(edge[id0], edge[idt]);
					boundary = max(boundary,1);
					insertArc(id0, idt, boundary);
				}
			}
		}
	}

	double cosineDist(unsigned short *A, unsigned short *B)
	{
		double mul=0, d_a=0, d_b=0;

		for(int i=0; i<m_chnl; ++i)
		{
			mul += A[i]*B[i];
			d_a += A[i]*A[i];
			d_b += B[i]*B[i];	
			//mexPrintf("A-%i B-%i mul-%i", A[i],B[i],mul);
			//mexEvalString("drawnow;");
		}
		// if ((unsigned short)sqrt((double)d_a*d_b)==0)
		// {
		// 	return 65334;
		// }

		// mexPrintf("mul-%g d_a-%g d_b-%g cos-%g", mul,d_a,d_b,(1-(mul)/(sqrt(d_a*d_b))));
		// mexEvalString("drawnow;");

		return (1-(mul)/(sqrt(d_a*d_b))); 
	}

		void growRegionJoin()
	{

		// // radix sort. Hace sorting de los arcs segun su dist!
		// radixSortLSD(m_dist,m_dtemp,m_minarc,m_arctmp,m_vexnum);
		// // if (m_vexnum<10)
		// // {

		// 	// mexPrintf("i_to_vex_Num = %i.\n",(m_vexnum*m_selPerc)/100);
		// 	// mexEvalString("drawnow;");
		// 	for (int i=0; i<ceil((m_vexnum*m_selPerc)/100); ++i)
		// 	{
				
		// 		ArcBox *arc = m_minarc[i];
				
		// 		int pu = findset(arc->u);
		// 		int pv = findset(arc->v);
		// 		// mexPrintf("dist = %i, pu = %i, arcpu = %i.\n", m_dist[i], pu, arc->u);
		// 		// mexEvalString("drawnow;");
		// 		if (pu == pv) continue;

		// 		m_treeu[m_treeSize]   = pu;
		// 		m_treev[m_treeSize++] = pv;
		// 		// m_countVex[pu]++;
		// 		// m_countVex[pv]++;
		// 		// mexPrintf("pu=%i, counting=%i \n",pu ,m_countVex[pu]);
		// 		// mexEvalString("drawnow;");
		// 		if (pu < pv) m_parent[pv] = pu;
		// 		else         m_parent[pu] = pv;
		// 	}
			
		// 	// update parent and label
		// 	m_regionnum = 0;
		// 	for (int i=0; i<m_vexnum; ++i)
		// 	{
		// 		int u = m_vex[i];
		// 		int pu = findset(u);
				
		// 		// mexPrintf("u = %i, pu = %i.\n", u, pu);
		// 		// mexEvalString("drawnow;");

		// 		if (pu == u)
		// 		{
		// 			m_label[u] = m_regionnum;
		// 			// create a new component
		// 			m_vex[m_regionnum] = u;
		// 			m_frtarc[m_regionnum] = NULL; //Eliminate arc
		// 			++m_regionnum;
		// 		}
		// 		else
		// 		{
		// 			m_label[u] = m_label[pu];
		// 			// update color
		// 			int s1 = m_size[u]; int s2 = m_size[pu]; int s = s1+s2;
		// 			m_size[pu] = s;

		// 			// Sparseness color of both. [defined by Hoyer]
		// 			int sp1sum=1; int sp1sum2=1; int sp2sum=1; int sp2sum2=1;
		// 			// mexPrintf("Sparseness Color start:\n");
		// 		    // mexEvalString("drawnow;");
		// 			for (int ii=0; ii<m_chnl; ++ii) 
		// 			{
		// 				sp1sum =  sp1sum+m_color[m_chnl*pu+ii];
		// 				sp1sum2 =  sp1sum+m_color[m_chnl*pu+ii]*m_color[m_chnl*pu+ii];
		// 				sp2sum =  sp1sum+m_color[m_chnl*u+ii];
		// 				sp2sum2 =  sp1sum+m_color[m_chnl*u+ii]*m_color[m_chnl*u+ii];						
		// 			}

		// 			// mexPrintf("sp1sum2 = %i, sp2sum2 = %i.\n", sp1sum2, sp2sum2);
		// 		    // mexEvalString("drawnow;");
		// 			sp1sum2 = sqrt(sp1sum2); sp2sum2 = sqrt(sp2sum2);					
		// 			if (((sqrt(m_chnl)-sp1sum/sp1sum2)/(sqrt(m_chnl)-1))<((sqrt(m_chnl)-sp2sum/sp2sum2)/(sqrt(m_chnl)-1)))
		// 			{
		// 				for (int ii=0; ii<m_chnl; ++ii) m_color[m_chnl*pu+ii] = m_color[m_chnl*u+ii];
		// 			}
		// 			//arc->confidence = s+arc->confidence;
		// 		}
				
		// 	}

		// 	m_vexnum = m_regionnum;
			
		// 	// mexPrintf("vex_num %i.\n", m_vexnum);
		// 	// mexEvalString("drawnow;");
	}

};

#endif