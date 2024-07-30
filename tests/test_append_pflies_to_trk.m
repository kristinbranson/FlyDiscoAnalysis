load('append_flies_to_trk_inputs.mat', 'ft_trk', 'trx') ;
input_trk = ft_trk ;

output_trk = append_pflies_to_trk(input_trk, trx)
