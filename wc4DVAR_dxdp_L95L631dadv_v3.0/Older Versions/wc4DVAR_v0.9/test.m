for i=1:n 
xl(:,i+1) = Mnl(xl(:,i),fnl,schnl,dt,1);    
[ddx(:,i+1),~,Mm] = Mtlmtx(ddx(:,i),xl(:,i),ftl,fnl,schtl,schnl,dt,1);
              M{1,i} = Mm;
end

[ddx(:,i+1),~,Mm] = Mtlmtx(ddx(:,i),xl(:,i),ftl,fnl,schtl,schnl,dt,2);