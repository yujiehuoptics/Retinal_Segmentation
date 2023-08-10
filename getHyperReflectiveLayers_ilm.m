function pathsf = getHyperReflectiveLayers_ilm(inputImg,constants)

%%%胡宇捷编写
%$Revision: 1.1 $ $Date: 2013/09/15 21:00$ $Author: Pangyu Teng $
    % parameter for smothing the images.     %图像平滑参数
%     params.filter0Params = [5 5 1];
    params.filterParams = [20 20 2];  

%for a very smooth image, a "broad stroke" of the image
inputImg = imfilter(inputImg,fspecial('gaussian',params.filterParams(1:2),params.filterParams(3)),'replicate');  

if nargin < 1
    display('requires at least 1 input (findHyperReflectiveZones.m)');
    return;
end

if nargin == 1
    %initiate parameters
    constants.shrinkScale = 0.2;     %缩放因子
    constants.offsets = -20:20;      %偏移量
end

isPlot = 0;

%shrink the image.
szImg = size(inputImg);
procImg = imresize(inputImg,constants.shrinkScale,'bilinear');   %图像缩小
% procImg =inputImg;
%create adjacency matrices
[adjMatrixW, adjMatrixMW, adjMX, adjMY, adjMW, adjMmW, newImg] = getAdjacencyMatrix(procImg);    %计算邻接矩阵

%create roi for getting shortestest path based on gradient-Y image.
[gx, gy] = gradient(newImg);     %newImg的水平、垂直梯度向量 newImg已经将图片左右各加两列了
szImgNew = size(newImg);
roiImg = zeros(szImgNew);
roiImg(gy > mean(gy(:))) =1 ;  %把其中垂直梯度值大于平均值的元素置为1，由黑到白明显的是1
figure();
imshow(roiImg,[]);
% find at least 2 layers    %最少找出两层
path{1} = 1;
count =1
while ~isempty(path) && count <= 2

    %add columns of one at both ends of images  %在原图像基础上增添第一列和最后一列，元素为1
    roiImg(:,1)=1;
    roiImg(:,end)=1;
    
    % include only region of interst in the adjacency matrix
    includeX = ismember(adjMX, find(roiImg(:) == 1));           %返回原图中值垂直梯度等1的索引值，且对应位置为1   find(roiImg(:) == 1找到索引值位置
    includeY = ismember(adjMY, find(roiImg(:) == 1));
    keepInd = includeX & includeY;
    
    % compile adjacency matrix
    adjMatrix = sparse(adjMX(keepInd),adjMY(keepInd),adjMmW(keepInd),numel(newImg(:)),numel(newImg(:)));
     %sparse是生成一个稀释矩阵，S=sparse(i,j,s,m,n)是用向量i，j和s生成一个m×n的稀释矩阵，i和j是行下标和列下标
    % S = sparse(i,j,s,m,n,nzmax)   S(i(k),j(k)) = s(k)
    
    % get layer going from dark to light   %这是采用封装好的Dijkstra 算法用于最短路径问题
    %sparse([1,2,3,4],[1,2,3,4],[0,0,1,1],5,5,6)  ans =(3,3) 1 (4,4) 1   
    %adjMX(keepInd)是22616
    %为8幅图片的有效数据，每幅图片为10710，所以为10710的稀疏矩阵，矩阵每一行的非零元素表示图像中的一个点与周围的坐标关系以及权重   
    [ dist,path{1} ] = graphshortestpath( adjMatrix, 1, numel(newImg(:)));
    
    if ~isempty(path{1})    %如果没有找到
                        
        % get rid of first few points and last few points
        [pathX,pathY] = ind2sub(szImgNew,path{1});           %pathX,pathY 都为行向量  

        pathX = pathX(gradient(pathY)~=0);        %只取pathY中垂直梯度分量不为0的pathX值和pathY值   这种情况多的话就会出现直线或者长的尖的斜线
        pathY = pathY(gradient(pathY)~=0);
        
        %block the obtained path and abit around it                 %在坐标区域重复pathX，pathY  行复制41次，列复制41次
        pathXArr = repmat(pathX,numel(constants.offsets));
        pathYArr = repmat(pathY,numel(constants.offsets));
        for i = 1:numel(constants.offsets)
            pathYArr(i,:) = pathYArr(i,:)+constants.offsets(i);    %将-20：20加入到对应的每行
        end
        
        pathXArr = pathXArr(pathYArr > 0 & pathYArr <= szImgNew(2));
        pathYArr = pathYArr(pathYArr > 0 & pathYArr <= szImgNew(2));
        
        pathArr = sub2ind(szImgNew,pathXArr,pathYArr);
        roiImg(pathArr) = 0;
        
        paths(count).pathX = pathX;
        paths(count).pathY = pathY;

        if isPlot;
            subplot(1,3,1);
            imagesc(inputImg);
            subplot(1,3,2);
            imagesc(gy);        
            subplot(1,3,3);
            imagesc(roiImg);
            drawnow;
            pause;
        end
        
    end % of ~empty
    count = count + 1;
end


if ~exist('paths','var')
    paths = {};
    keyboard;
    return;
end % if exist

%format paths back to original size
for i = 1:numel(paths)    
    [paths(i).path, paths(i).pathY, paths(i).pathX] = resizePath(szImg, szImgNew, constants, paths(i).pathY, paths(i).pathX);    
    paths(i).pathXmean = nanmean(paths(i).pathX);
    paths(i).name = [];
end
    

%name each path (numel(paths) should equal to 2)
% if numel(paths) ~= 2
%     paths = {};
%     display('error');
%     return;
% end

% based on the mean location detemine the layer type.
if paths(1).pathXmean < paths(2).pathXmean                %根据层位置得到层名字
    paths(1).name = 'ilm';
    paths(2).name = 'isos';
    pathsf = paths(1);
else
    paths(1).name = 'isos';    
    paths(2).name = 'ilm';    
    pathsf = paths(2);
end
% 
% 
% if isPlot;
%     imagesc(inputImg);
%     axis image; colormap('gray');
%     hold on;
%     for i = 1:numel(paths)
%         cola = rand(1,3);
%         plot(paths(i).pathY,paths(i).pathX,'r-','linewidth',3);
%         text(paths(i).pathY(end),paths(i).pathX(end)-15,paths(i).name,'color',rand(1,3));
%         drawnow;
%     end
%     
%     hold off;
end