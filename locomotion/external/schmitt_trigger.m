%SCHMITT TRIGGER
function [y] = schmitt_trigger(x,tL,tH)

% x=rand(1,100);
% x=conv(x,ones(10,1)/10);

if tL > tH
    disp('error: low threshold is higher then high threshold')
    return
end


 limit=0;
if nargin<3
   disp('number of inputs arguments should be three');
end
   N=length(x);
  
   y=[length(x)];
  

   
   for i=1:N
       
       
      if ( limit ==0)
          
          y(i)=0;
          
      elseif (limit == 1)
           
           y(i)=1;
           
      end
       
       
      if (x(i)<=tL)
          limit=0; 
            y(i)=0;
          
      elseif( x(i)>= tH)         
          limit=1;  
          y(i)=1;
              
      end
      
         
   end
  
%   plot(x,'r','DisplayName','plot of x','LineWidth',1.5); hold on;
%   plot(y,'blue','DisplayName','plot of y','LineWidth',3); hold off;
%   legend('show');

% modified Alice Robie 2023
%   Copyright (c) 2013, jamali omary
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
