function tri = surf2tri ( p, t )                
                %
                %p, t, pause 
                faces = [  t(:,[1,2,3]);
                           t(:,[1,2,4]);
                           t(:,[1,3,4]);
                           t(:,[2,3,4])];
                node4 = [t(:,4);t(:,3);t(:,2);t(:,1)];
                faces = sort(faces,2);
                [~,iindx,jindx] = unique(faces,'rows');
                vec = histc(jindx,1:max(jindx)); 
                qx = find(vec==1);
                tri = faces(iindx(qx),:);
                node4 = node4(iindx(qx));
                %
                % Orientation
                %
                v1 = p(tri(:,2),:)-p(tri(:,1),:);
                v2 = p(tri(:,3),:)-p(tri(:,1),:);
                v3 = p(node4,:)-p(tri(:,1),:);
                iindx = find(dot(cross(v1,v2,2),v3,2)>0);
                tri(iindx,[2,3]) = tri(iindx,[3,2]);
   end
