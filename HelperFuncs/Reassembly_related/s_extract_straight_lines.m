function [ lines, lines_numel ] = s_extract_straight_lines( poly, thresh )
%S_EXTRACT_STRAIGHT_LINES find the largest straight lines in a given
%polygon where the error is not largest than thresh
% input: poly: 2xN
EPSi=1e-13;
N=size(poly,2);
if size(poly,1)~=2
    error('wrong input. poly size=[%d,%d]',size(poly,1),size(poly,2));
end
lines=zeros(0,2);
lines_numel=zeros(0,1);
V_visited=false(1,N);
INITIAL_STEP=10;

si=1;
ei=cyc_next(si,N,INITIAL_STEP);

direction = 1;
while true %will break when si is empty
%     V_visited(si)=true;
%     V_visited(ei)=true;
    if si<ei
        x=poly(1,si:ei);
        y=poly(2,si:ei);
        V_visited(si:ei)=true;
    else
        x=[poly(1,si:N),poly(1,1:ei)];
        y=[poly(2,si:N),poly(2,1:ei)];
        V_visited(1:ei)=true;
        V_visited(si:N)=true;
    end
    [~,errorsq,~,xfit,yfit]=points2linecoeff(x,y);
    
    
    
    if (direction>0 && errorsq>thresh) ||...
       (direction>0 && V_visited(cyc_next(ei,N)))
%    figure(1);clf(1);
%     plot(poly(1,:),poly(2,:),'r',xfit,yfit,'b');
%     legend('poly','current line');
%     title(sprintf('error: %0.4f, npoints: %d',errorsq, numel(xfit)));
        if errorsq>thresh
            ei=cyc_prev(ei,N);
        end
        direction=-1;
    elseif (direction<0 && errorsq>thresh) ||...
       (direction<0 && V_visited(cyc_prev(si,N)))
%    figure(1);clf(1);
%     plot(poly(1,:),poly(2,:),'r',xfit,yfit,'b');
%     legend('poly','current line');
%     title(sprintf('error: %0.4f, npoints: %d',errorsq, numel(xfit)));
        if errorsq>thresh
            si=cyc_next(si,N);
        end
        currline_cyc_numel = cyc_numel(si,ei,N);
        if currline_cyc_numel>0.1*N
            lines(end+1,:)=[si,ei];
            lines_numel(end+1)=currline_cyc_numel;
        end
        direction = 1;
        % find the next first vertex that wasnt visited and the following vertex also not visited
        si=find(diff([V_visited,V_visited(1:INITIAL_STEP)],INITIAL_STEP)==0 & ~V_visited,1);
        if isempty(si)
            break;
        end
        ei=cyc_next(si,N,INITIAL_STEP-1);
    end
    if direction > 0
        ei=cyc_next(ei,N);
    else
        si=cyc_prev(si,N);
    end
end
end
function c=cyc_prev(a,N,step)
if ~exist('step','var')
    step=1;
end
c=a-step;
if c<1
    c=c+N;
end
end
function c=cyc_next(a,N,step)
if ~exist('step','var')
    step=1;
end
c=a+step;
if c>N
    c=c-N;
end
end
function c=cyc_numel(a,b,N)
if a<b
    c=b-a+1;
else
    c=N-a+1+b;
end
end
