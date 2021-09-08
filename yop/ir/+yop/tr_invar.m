function tri = tr_invar(nbc)
% TR_INVAR - Transcription invariance transformation
% Separates expressions that are invariant in the parametrization to 
% changes in the independent variable during the transcription, from those
% that needs a different parametrization depending on the value of the
% independent variable. It makes the transcription easier as it suffices to
% test 'all(is_transcription_invariant)' and make a binary decision.

tri = yop.tri_data();

for k=1:length(nbc.eq)
    r = nbc.eq{k};
    invariant = is_transcription_invariant(r.lhs);
    if all(invariant) || all(~invariant)
        tri.add_eq(r);  
    else
        tri.add_eq(yop.get_subrelation(r, invariant));
        tri.add_eq(yop.get_subrelation(r, ~invariant));
    end
end

for k=1:length(nbc.ieq)
    r = nbc.ieq{k};
    invariant = is_transcription_invariant(r.lhs);
    if all(invariant) || all(~invariant)
        tri.add_ieq(r);
    else
        tri.add_ieq(yop.get_subrelation(r, invariant));
        tri.add_ieq(yop.get_subrelation(r, ~invariant));
    end
end

end