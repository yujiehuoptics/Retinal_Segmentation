%     {{Caserel}}
%     Copyright (C) {{2013}}  {{Pangyu Teng}}
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%


function paths = getHyperReflectiveLayers(inputImg,constants)
%$Revision: 1.1 $ $Date: 2013/09/15 21:00$ $Author: Pangyu Teng $

if nargin < 1
    display('requires at least 1 input (findHyperReflectiveZones.m)');
    return;
end

%%%%%%%%%�ñ��˼�벻��%%%%%%%%%%%%%%%%%%
if nargin == 1
    %initiate parameters
    constants.shrinkScale = 0.2;
    constants.offsets = -20:20;        %Ϊʲô��-20��20,����20������
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isPlot = 1;

%shrink the image.                  %ͼ����Ϊԭ����1/5
szImg = size(inputImg);
procImg = imresize(inputImg,constants.shrinkScale,'bilinear');

%create adjacency matrices
[adjMatrixW, adjMatrixMW, adjMX, adjMY, adjMW, adjMmW, newImg] = getAdjacencyMatrix(procImg);

%create roi for getting shortestest path based on gradient-Yimage.���ڴ�ֱ�ݶ�ͼ��ΪѰ�����·������roi����
[gx, gy] = gradient(newImg);              %newImg��ˮƽ����ֱ�ݶ����� newImg�Ѿ���ͼƬ���Ҹ���������
szImgNew = size(newImg);
roiImg = zeros(szImgNew);
roiImg(gy > mean(gy(:))) =1 ;             %�����д�ֱ�ݶ�ֵ����ƽ��ֵ��Ԫ����Ϊ1���ɺڵ������Ե���1

% find at least 2 layers                %�����ҳ�����
path{1} = 1;
count = 1;
while ~isempty(path) && count <= 2

    %add columns of one at both ends of images      %��ԭͼ�������������һ�к����һ�У�Ԫ��Ϊ1
    roiImg(:,1)=1;
    roiImg(:,end)=1;
    
    % include only region of interst in the adjacency matrix
    includeX = ismember(adjMX, find(roiImg(:) == 1));           %����ԭͼ��ֵ��ֱ�ݶȴ���1������ֵ���Ҷ�Ӧλ��Ϊ1   find(roiImg(:) == 1�ҵ�����ֵλ��
    includeY = ismember(adjMY, find(roiImg(:) == 1));          %���ر任��ͼ�������ֵ���Ҷ�Ӧλ��Ϊ1
    keepInd = includeX & includeY; %Ϊʲôԭͼ�ͱ任���ͼ����Ҫ�����ݶ�ֵ����1 ��������Ԫ�ص��ݶ�ֵ����Ҫ��������
    
    % compile adjacency matrix
    adjMatrix = sparse(adjMX(keepInd),adjMY(keepInd),adjMmW(keepInd),numel(newImg(:)),numel(newImg(:)));  
    %sparse������һ��ϡ�;���S=sparse(i,j,s,m,n)��������i��j��s����һ��m��n��ϡ�;���i��j�����±�����±�
    % S = sparse(i,j,s,m,n,nzmax)   S(i(k),j(k)) = s(k)
    % get layer going from dark to light   %���ǲ��÷�װ�õ�Dijkstra �㷨�������·������
    %sparse([1,2,3,4],[1,2,3,4],[0,0,1,1],5,5,6)  ans =(3,3) 1 (4,4) 1   
    %adjMX(keepInd)��22616
    %Ϊ8��ͼƬ����Ч���ݣ�ÿ��ͼƬΪ10710������Ϊ10710��ϡ����󣬾���ÿһ�еķ���Ԫ�ر�ʾͼ���е�һ��������Χ�������ϵ�Լ�Ȩ��
    [ dist,path{1} ] = graphshortestpath( adjMatrix, 1, numel(newImg(:)));  %��һ�ڵ㵽�����ľ��룬Ϊб�� path(1)Ϊ���������ֵ
    
    if ~isempty(path{1})
                        
        % get rid of first few points and last few points
        [pathX,pathY] = ind2sub(szImgNew,path{1});       %pathX,pathY ��Ϊ������  

        pathX = pathX(gradient(pathY)~=0);          %ֻȡpathY�д�ֱ�ݶȷ�����Ϊ0��pathXֵ��pathYֵ   ���������Ļ��ͻ����ֱ�߻��߳��ļ��б��
        pathY = pathY(gradient(pathY)~=0);
        
        %block the obtained path and abit around it       %�����������ظ�pathX��pathY  �и���41�Σ��и���41��           
        pathXArr = repmat(pathX,numel(constants.offsets));  %
        pathYArr = repmat(pathY,numel(constants.offsets));
        for i = 1:numel(constants.offsets)
            pathYArr(i,:) = pathYArr(i,:)+constants.offsets(i);   %��-20��20���뵽��Ӧ��ÿ��
        end
        
        pathXArr = pathXArr(pathYArr > 0 & pathYArr <= szImgNew(2));
        pathYArr = pathYArr(pathYArr > 0 & pathYArr <= szImgNew(2));
        
        pathArr = sub2ind(szImgNew,pathXArr,pathYArr);      
        roiImg(pathArr) = 0;        %roilImgΪ105*102��pathArrΪ154242*1    %ȥ��·���ϵ�����20�����ص�
        
        paths(count).pathX = pathX;
        paths(count).pathY = pathY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ҫ%%%%%%%%%%%%%%%
        if isPlot;
            subplot(1,3,1);
            imagesc(inputImg);
            subplot(1,3,2);
            imagesc(gy);        
            subplot(1,3,3);
            imagesc(roiImg);
            drawnow;
%             pause;
        end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    end % of ~empty
    count = count + 1;
end
paths;

%%%%%%%%%%%��paths������ʱ�����Գ�����
if ~exist('paths','var')
    paths = {};
    keyboard;
    return;
end % if exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%format paths back to original size   %��ʽ��·��Ϊ����ĳߴ�
for i = 1:numel(paths)    
    [paths(i).path, paths(i).pathY, paths(i).pathX] = resizePath(szImg, szImgNew, constants, paths(i).pathY, paths(i).pathX);    
    paths(i).pathXmean = nanmean(paths(i).pathX);
    paths(i).name = [];
    
end

%name each path (numel(paths) should equal to 2��)   
if numel(paths) ~= 2
    paths = {};
    display('error');
    return;
end

%based on the mean location detemine the layer type.   %����λ��ƽ��ֵ�ֳ������
if paths(1).pathXmean < paths(2).pathXmean
    paths(1).name = 'ilm';
    paths(2).name = 'isos';
else
    paths(1).name = 'isos';    
    paths(2).name = 'ilm';    
end


if isPlot;                                            %���ֲ��� 
    imagesc(inputImg);
    axis image; colormap('gray');
    hold on;
    for i = 1:numel(paths)
%         cola = rand(1,3);
        plot(paths(i).pathY,paths(i).pathX,'r-','linewidth',3);
        text(paths(i).pathY(end),paths(i).pathX(end)-15,paths(i).name,'color',rand(1,3));
        drawnow;
    end
    hold off;
end