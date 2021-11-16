%% Lacunarity Index Computing With Gliding Box Algorithm
% Author : Amirhossein Samavi(amirhosseinsamavi79@gmail.com)
function La = LACUNARITY(Ipath)
%% prepare image
I=imread("14200110.bmp");
I=im2bw(I,0.4);
I=-I+1;
%I = round(rand(8,16));
%I = [1,1,0,1,1,1,0,1,0,1,1,0;0,0,0,0,0,1,0,0,0,1,1,1;0,1,0,1,1,1,1,1,0,1,1,0;1,0,1,1,1,0,0,0,0,0,0,0;1,1,0,1,0,1,0,0,1,1,0,0;0,1,0,1,1,0,0,1,0,0,1,0;0,0,0,0,0,1,1,1,1,1,1,1;0,1,1,0,0,0,1,1,1,1,0,0;0,1,1,1,0,1,1,0,1,0,0,1;0,1,0,0,0,0,0,0,0,1,1,1;0,1,0,1,1,1,0,1,1,0,1,0;0,1,0,0,0,1,0,1,1,1,0,1]

%% Define box sizes and other variables
[rows,cols] = size(I);
boxSizes = [1,2,4,8,16,32,64];
maxSize = max(boxSizes);
len = length(boxSizes);
% N(r) : total number of boxes with r size and M mass
N = zeros(1,len);
Q = zeros((maxSize^2)+1,len);
Z1 = zeros((maxSize^2)+1,len);
Z2 = zeros((maxSize^2)+1,len);
M = transpose(0:maxSize^2);
squareM = M.^2;
for k=1:len
    currentSize = boxSizes(k)
    bnx = rows-currentSize+1;
    bny = cols-currentSize+1;
    N(k) = bnx*bny;
    for i=1:bnx
        for j=1:bny
            mass = sum(sum(I(i:i+currentSize-1, j:j+currentSize-1)));
            % mass+1 expression because of Q(0) is undefined
            Q(mass+1 ,k) = Q(mass+1 ,k) + 1;
        end
    end
    Q(:,k) = Q(:,k)./N(k);
    %check the result sum of Q always must be 1
    checkAllIsOne = sum(Q)
    %Z(n) = sigma(M^n * Q(M,r))
    Z1(:,k) = Q(:,k).*M;
    Z2(:,k) = Q(:,k).*squareM;
end
Lacunarity = sum(Z2)./sum(Z1).^2
plot(boxSizes,Lacunarity,'r-x');hold on
xlabel('box sizes','FontSize',14);
ylabel('Lacunarity index','FontSize',14);
