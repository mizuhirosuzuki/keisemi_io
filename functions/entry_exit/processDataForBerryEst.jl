function processDataForBerryEst(df)
    
    data_processed = df[:, :];
    data_processed[!, :TieEntryOrder] = rand(nrow(df));
    
    transform!(
        groupby(data_processed, :CityCode),
        nrow => :numPotenHos,
        :MRIOwnDum => sum => :numEntryObs
    );
    sort!(
        data_processed, 
        [:CityCode, :LogNumBeds, :TieEntryOrder], 
        rev = [false, true, true]
        );
    transform!(groupby(data_processed, :CityCode), :numPotenHos => (x -> 1:length(x)) => :EntryOrderId);
    
    return data_processed
        
end
