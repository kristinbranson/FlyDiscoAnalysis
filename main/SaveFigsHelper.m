function SaveFigsHelper(hfig,filename,dirname)


ISDIR = nargin >= 3 && ~isempty(dirname);
if ISDIR,
  figfilename = fullfile(dirname,[filename,'.fig']);
else
  figfilename = [filename,'.fig'];
end  

saveas(hfig,figfilename,'fig');
savefig([filename,'.pdf'],hfig,'pdf');
%saveas(hfig,[filename,'.pdf'],'pdf');
if ISDIR,
  unix(sprintf('mv %s.pdf %s/%s.pdf',filename,dirname,filename));
end
saveas(hfig,[filename,'.jpg'],'jpg');
%savefig([filename,'.jpeg'],hfig,'jpeg');
if ISDIR,
  unix(sprintf('mv %s.jpg %s/%s.jpeg',filename,dirname,filename));
end