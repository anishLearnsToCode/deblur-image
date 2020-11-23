function [s,tem]=kurtosis_p_k(B)


    B = double(B);
%     g=B;
% g = fspecial('unsharp');
% B = conv2(B,g,'same');
% g = fspecial('log');
% B = conv2(B,g,'same');    

p_size=80;

[r,c]=size(B);
padB=ones(r+p_size,c+p_size);
padB((p_size/2+1):end-(p_size/2),(p_size/2+1):end-(p_size/2))=B;
sz=size(((p_size/2+1):2:r+(p_size/2)),2)*size(((p_size/2+1):2:c+(p_size/2)),2);
% temp=cell(sz);
% tem=cell(sz);
pad=B;
in=1;
% 'finding kurt'
for i=1:20:r-(p_size)
    i; 
    for j=1:20:c-(p_size)
        temp=padB(i:i+(p_size),j:j+(p_size));
        t(in,:)=temp(:);
        k(in)=kurtosis(temp(:));
                in=in+1;

    end
end

mini=min(min(k));
[s,ind]=sort(k(:),'descend');
size(s);
size(t);
% 'patch arrangement'
for i=1:size(s)
     tem(i,:)=t(ind(i),:);
end
% 
% % 
% % [x,y]=find(k==s(10) );
% % if x<(p_size/2)
% %     x=(p_size/2);
% % elseif x>r-(p_size/2)
% %     x=r-(p_size/2);
% % end
% % if y<(p_size/2)
% %     y=(p_size/2);
% % elseif y>c-(p_size/2)
% %     y=c-(p_size/2);
% % end
% % 
% % p=B(x-(p_size/2-1):x+(p_size/2),y-(p_size/2-1):y+(p_size/2));
% 
% % [U,S,V]=svd(B);
% % S=diag(S);
% % [S,ind]=sort(S(:),'descend');
% % size(s)
% % size(t)
% % for i=1:size(s)
% %      U1(i,:)=U(ind(i),:);
% % end
% 
% for i=1:size(s)
%     ke=reshape(t(i,:),81,81);
% %     ke=power((ke),-4);
%     co=conv2(B,ke/max(ke(:)));
%     n(i) = norm(co);
% end
% figure()
% plot(s/max(s))
% hold on
% plot((n/max(n)))
% hold on
end