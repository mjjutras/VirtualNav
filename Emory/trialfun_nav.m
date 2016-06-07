function trl = trialfunclrchng(cfg)

% read the header
hdr = ft_read_header(cfg.dataset);
% read the events
event = ft_read_event(cfg.dataset);

numevt = length(event);
for evtlop = 1:numevt
  mrk.val(evtlop) = event(evtlop).value;
  mrk.tim(evtlop) = event(evtlop).sample;
end


trl= [];
trldum = [];
for rptlop=1:length(mrk.val)
    if mrk.val(rptlop)==100
        if isempty(trldum)
            trldum(1) = mrk.tim(rptlop);
        else
            trldum(2) = mrk.tim(rptlop)-1;            
            trl = [trl; [[trldum] 0]];
            trldum(1) = mrk.tim(rptlop);
        end
    end
end



           
