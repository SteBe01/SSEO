function plotObj(r,offset,color)
    if nargin==2
        color=[0.4660 0.6740 0.1880];
    end
    % Make unit sphere
    [x,y,z] = sphere;
    % Scale to desire radius.
    x = x * r;
    y = y * r;
    z = z * r;
    % Plot as surface translated in new centre.
    surf(x+offset(1),y+offset(2),z+offset(3),'FaceColor',color,'EdgeColor','none'); 
    axis equal;
    hold on;
    
end

