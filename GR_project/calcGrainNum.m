function [numOfGrains,numOfBud] = calcGrainNum(imagePath)
I=imread(imagePath);
% imtool(I);
I_tmp = I;
I_tmp1 = I;
I_tmp2 = I;
I=double(I);
I_gray=I(:,:,1);
I_gray1=I(:,:,2);
size1=size(I);
I_column=reshape(I,size1(1)*size1(2),3);
thr1=60;thr2=60;
I1=I_column(:,1)-I_column(:,3);I1(I1>=thr1)=30;
I2=I_column(:,2)-I_column(:,3);I2(I2>=thr2)=30;
I4=I_column(:,1)+I_column(:,2)+I_column(:,3);
I1=[I1,I2];
matrix1=[40 40;0 0]; %定义2个初始聚类中心
idx=kmeans(I1,2,'Start',matrix1); %利用Start指定初始聚类中心
idx=idx-1;
temp1=mean(I4(idx==0));
temp2=mean(I4(idx==1));
if temp1>temp2
    idx=~idx;
end
idx_reshape=reshape(idx,size1(1),size1(2));


I1=idx_reshape;
I1 = bwareaopen(I1,80,4);%删除较小面积的连通图
I1=imfill(I1,'holes');

% [labeled,numObjects]=bwlabel(I1,4); %确定图像中的连通区域个数
labeled = bwlabeln(I1,4);
numObjects=max(max(labeled));

%缩小每个连通区域的谷粒
each_grain_area_histogram=zeros(numObjects,256);   %统计连通区域红色分量直方图
for i=1:size1(1)
    for j=1:size1(2)
        temp=labeled(i,j);
        if temp~=0
            temp1=I_gray(i,j);
            temp1=temp1+1;
            each_grain_area_histogram(temp,temp1)=each_grain_area_histogram(temp,temp1)+1;            
        end
    end
end
each_grain_area_histogram_inflection_position=zeros(numObjects,1);
x=1:256;
d=5;
initial1=150;
for i=1:numObjects %对直方图进行平滑
     h=fspecial('gaussian',[1,100],7.5);
    each_grain_area_histogram(i,:)=imfilter(each_grain_area_histogram(i,:),h);
    %确定直方图的趋势线
    y=each_grain_area_histogram(i,:);
    r=polyfit(x,y,d);
    yv=r(1)*x.^5+r(2)*x.^4+r(3)*x.^3+r(4)*x.^2+r(5)*x+r(6);%直方图的趋势线
    [max1,indx1]=max(yv(initial1:255));
    indx1=indx1+initial1-2;
    for k=indx1:-1:2
        temp1=yv(k);%沿最大值向左下角找第一个谷点
        temp2=yv(k-1);
        temp3=yv(k+1);        
        if temp1<=temp2 && temp1<=temp3 
            break;
        end
    end
    each_grain_area_histogram_inflection_position(i,1)=k;
end
for i=1:size1(1)
    for j=1:size1(2)
        temp=labeled(i,j);
        if temp~=0
            if each_grain_area_histogram_inflection_position(temp,1)~=0 ...
                && I_gray(i,j)<each_grain_area_histogram_inflection_position(temp,1)
                I1(i,j)=0;
                labeled(i,j)=0;
            end                
        end
    end
end
 I1 = bwareaopen(I1,100,4);%删除较小面积的连通图

[labeled,numObjects]=bwlabel(I1,4); %确定图像中的连通区域个数

each_grain_area=zeros(numObjects,1); %统计每个谷粒连通区域的面积
for i=1:size1(1)
    for j=1:size1(2)
        temp1=labeled(i,j);
        if temp1~=0
            each_grain_area(temp1)=each_grain_area(temp1)+1;
        end
    end
end
[sort_each_grain_area indx]=sort(each_grain_area);

h1=fspecial('gaussian',3,0.65);
smooth_sort_each_grain_area=imfilter(sort_each_grain_area,h1); %面积分布曲线图高斯平滑

sort_each_grain_area_grad=gradient(smooth_sort_each_grain_area,1);

sort_each_grain_area_grad2=gradient(sort_each_grain_area_grad,1);%梯度曲线再求梯度

start0=0; end0=0;
for i=2:numObjects-1
    if sort_each_grain_area_grad2(i-1)<0&&sort_each_grain_area_grad2(i+1)>0
        start0=i;
        break;
    end
end
for i=start0+1:numObjects-3
    temp=(sort_each_grain_area(i+3)-sort_each_grain_area(i))/sort_each_grain_area(i);
    if temp>0.2
        end0=i;
        break;
    end
end

if end0==0     %图像中的谷粒都是单个独立摆放的，无粘连
    end0=numObjects;
end

% h2=fspecial('gaussian',7,10);
% smooth_sort_each_grain_area_grad=imfilter(sort_each_grain_area_grad,h2); %梯度曲线图高斯平滑
% figure;plot(x,smooth_sort_each_grain_area_grad);
%确定最优单个谷粒面积
min_error=1000000;
for i=start0:end0
    total_error=0;
    total_number=0;
    for j=start0:end0
        temp1=fix(sort_each_grain_area(j)/sort_each_grain_area(i));%取整数部分
        temp3=sort_each_grain_area(j)/sort_each_grain_area(i);
        if (temp1==0 || temp1==1)
            total_error=total_error+abs((1-temp3));
            total_number=total_number+1;
        end
    end
    total_error=total_error/total_number;
    if total_error<min_error
        min_error=total_error;
        one_grain_area=sort_each_grain_area(i);
        opt_grain_index=indx(i);
    end
end
%谷粒计数
threshold1=0.5;
threshold2=0.4;
total1=0;
contain=zeros(numObjects,1);   %记录每个谷粒连通区域包含的谷粒数
contain1=zeros(numObjects,1);   %记录每个谷粒连通区域包含的谷粒数

for i=1:numObjects
    temp1=fix(each_grain_area(i)/one_grain_area);%取整数部分
    temp2=each_grain_area(i)/one_grain_area-temp1;%取小数部分
    contain1(i)=each_grain_area(i)/one_grain_area;
    if temp1==0&&temp2>threshold1
        total1=total1+1;
        contain(i,1)=1;
    end
    if temp1>0&&temp2>threshold1
        total1=total1+temp1+1;
        contain(i,1)=temp1+1;
    end
    if temp1>0&&temp2<threshold1
        total1=total1+temp1;
        contain(i,1)=temp1;
    end
end


optGrain=zeros(size1(1),size1(2));  
for i=1:size1(1)
    for j=1:size1(2)
        ttt=labeled(i,j);
        if ttt==opt_grain_index
%             g=I_gray(i,j);
            optGrain(i,j)=1;
        end
    end
end
%------------------统计萌发种子数目-----------------------------------------
white_threshold=160;%定义白色的阈值
white_connected=zeros(size1(1),size1(2));
for i=1:size1(1)
    for j=1:size1(2)
        if labeled(i,j)==0 && I_gray(i,j)>white_threshold&& I_gray1(i,j)>white_threshold
            white_connected(i,j)=1;
        else
            white_connected(i,j)=0;
        end
    end
end

% figure;imshow(white_connected)
min_bud_area=fix(one_grain_area/50);%确定最小胚芽面积
white_connected = bwareaopen(white_connected,min_bud_area,4);%删除较小面积的连通图
max_bud_area=fix(one_grain_area/3);%确定最大胚芽面积
white_connected = bwareaopen1(white_connected,max_bud_area,4);%删除较大面积的连通图

[white_labeled,white_numObjects]=bwlabel(white_connected,4); %确定图像中的连通区域个数
isWhiteConnectGrain=zeros(white_numObjects,1);   %标记白色连通区域是否与谷粒连通区域相连
for i=2:size1(1)-1
    for j=2:size1(2)-1
        temp1=white_labeled(i,j);
        if temp1==0
            continue;
        end
        if isWhiteConnectGrain(temp1)==1
            continue;
        end
        if temp1>0 && (labeled(i-1,j)>0 || labeled(i+1,j)>0 ...
            || labeled(i,j-1)>0 || labeled(i,j+1)>0)
            isWhiteConnectGrain(temp1)=1;   
        end
    end
end
new_white_connected=zeros(size1(1),size1(2)); %与谷粒相连的白色区域
for i=2:size1(1)-1
    for j=2:size1(2)-1
        temp1=white_labeled(i,j);
        if temp1==0
            continue;
        end
       if isWhiteConnectGrain(temp1)==1
            new_white_connected(i,j)=1;
       end
    end
end
new_white_labeled = bwlabeln(new_white_connected,8);
new_white_numObjects=max(max(new_white_labeled));
new_white_connected_perimeter=regionprops(new_white_labeled, 'Perimeter');   %白色连通区域周长

mark_touch_area=zeros(new_white_numObjects,1);    %标记每个胚芽所属的谷粒连通区域
touchPixelsNum=zeros(new_white_numObjects,1);     %标记每个白色连通区域与谷粒粘连处的像素点数目
isTouchPixels=zeros(size1(1),size1(2));         %标记白色像素点是否点与谷粒粘连

for i=2:size1(1)-1
    for j=2:size1(2)-1
        temp1=new_white_labeled(i,j);
        
        if temp1==0
            continue;
        end
        if temp1>0 && (labeled(i-1,j)>0 || labeled(i+1,j)>0 ...
            || labeled(i,j-1)>0 || labeled(i,j+1)>0)

            if labeled(i-1,j)>0
                mark_touch_area(temp1)=labeled(i-1,j);
                touchPixelsNum(temp1)=touchPixelsNum(temp1)+1;  %如果一个白色连通区域与多个不同的谷粒区域相连，算法只标记一个相连的谷粒连通区域
                isTouchPixels(i,j)=temp1;
                continue;
            elseif labeled(i+1,j)>0
                mark_touch_area(temp1)=labeled(i+1,j);
                touchPixelsNum(temp1)=touchPixelsNum(temp1)+1;
                isTouchPixels(i,j)=temp1;
                continue;
            elseif labeled(i,j-1)>0
                mark_touch_area(temp1)=labeled(i,j-1);
                touchPixelsNum(temp1)=touchPixelsNum(temp1)+1;
                isTouchPixels(i,j)=temp1;
                continue;
            else
                mark_touch_area(temp1)=labeled(i,j+1);
                touchPixelsNum(temp1)=touchPixelsNum(temp1)+1;
                isTouchPixels(i,j)=temp1;
            end    
        end
    end
end
 
div1=zeros(new_white_numObjects,1);  

mark_whiteRegion_isValid=zeros(new_white_numObjects,1);
for i=1:new_white_numObjects
       div1(i)=touchPixelsNum(i)/new_white_connected_perimeter(i).Perimeter;
      if div1(i)<0.4
        mark_whiteRegion_isValid(i)=1;
      end
end
new_white_connected3=zeros(size1(1),size1(2)); 
for i=1:size1(1)
    for j=1:size1(2)
        l=new_white_labeled(i,j);
        if l==0
            continue;
        end
        if(mark_whiteRegion_isValid(l)==0)
            new_white_connected3(i,j)=0;
        end
        if(mark_whiteRegion_isValid(l)==1)
            new_white_connected3(i,j)=1;
        end
    end
end
numOfGrains=total1;
numOfBud=sum(mark_whiteRegion_isValid);
if numOfBud>numOfGrains
    numOfBud=numOfGrains;
end
%谷粒彩色图像
for i=1:size1(1)
    for j=1:size1(2)
        tt=labeled(i,j);
        if tt>0
            I_tmp(i,j,1)=255;
            I_tmp(i,j,2)=255;
            I_tmp(i,j,3)=0;
        end
%         if tt==10
%             I_tmp(i,j,1)=255;
%             I_tmp(i,j,2)=0;
%             I_tmp(i,j,3)=0;
%         end
    end
end