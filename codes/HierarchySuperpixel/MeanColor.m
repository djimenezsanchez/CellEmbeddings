function [color, SuperPheno] = MeanColor(I, S)

% check image sizes
[hS, wS, cS] = size(S);
[hI, wI, cI] = size(I);
if hS~=hI || wS~= wI || cS~=1
    fprintf('Size error\n');
    return;
end
% segment index should start with 1
if min(S(:)) == 0
    S=S+1;
end

nLabel = max(S(:));
M = zeros(nLabel, cI);
A = zeros(nLabel, 1);
PositionsXY = zeros(nLabel,2);
for y=1:hI
    for x=1:wI
        id = S(y,x);
        for c=1:cI
            M(id,c) = M(id,c)+I(y,x,c);            
        end
        A(id) = A(id) + 1;
        PositionsXY(id,:) =PositionsXY(id,:)+[x,y]; 
    end
end

for id=1:nLabel
    PositionsXY(id,:) = round(PositionsXY(id,:)./A(id));
    for c=1:cI
        M(id,c) = M(id,c) / A(id);
    end
end

color = zeros(hI,wI,cI);
for y=1:hI
    for x=1:wI
        id = S(y,x);
        for c=1:cI
            color(y,x,c) = M(id,c);
        end
        % testingg
%         if x==PositionsXY(id,1) & y==PositionsXY(id,2)
%             color(y,x,3) = id; 
%         end
    end
end
SuperPheno=[PositionsXY,A,M/(-1+2^16)];
color = color;