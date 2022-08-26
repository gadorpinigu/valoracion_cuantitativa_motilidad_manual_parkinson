% % % SCRIPT ORIGINAL: PROYECTO TTKT
% ADAPTADO PARA TEST PARKINSON: GADOR PIÑEYRO, UNAV

function vid = TTK_videoinput(n)

k = imaqhwinfo('tisimaq_r2013');
Info = k.DeviceInfo.SupportedFormats;
vid = videoinput('tisimaq_r2013',n,Info{1,127});
%Configuramos las propiedades del trigger
triggerconfig(vid,'manual') 
%Una vez el vídeo empieza, el trigger se dispara continuamente hasta que se
%le da la orden de STOP
set(vid, 'TriggerRepeat', Inf);
%El valor default de FramesPerTrigger es 10
set(vid, 'FramesPerTrigger', 1);
 
end
