function Forces= ReadingForces(path, time)
    filename = path+string(time)+'\force' ;
    filetext = fileread(filename);

    expr  = 'Thermowell'; % Read from
    expr2 = 'fixed_Bottom'; % Read to
    begining = regexp(filetext,expr);
    ending = regexp(filetext,expr2);

    text = filetext(begining:ending);
    expr3= '(-?\d*([.]\d*)? -?\d*([.]\d*)? -?\d*([.]\d*)?)';
    text = regexp(text,expr3,'Match'); %Look for pattern of three numbers
    %Format correctly and string to double
    text= strrep(text,'  ','');
    text= text(~cellfun('isempty',text));
    text= split(text);
    text= squeeze(text);
    Forces= str2double(text);
end