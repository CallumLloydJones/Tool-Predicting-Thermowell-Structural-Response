function writeTractionFile(path, Traction, t)
    for ti=t
        mkdir(append(path,num2str(ti)));
        file = fopen( append(path, num2str(ti), '\traction'), 'wt' );
        mystring = '\';
        fprintf(file,'/*--------------------------------*- C++ -*----------------------------------*%s\n',mystring);
        fprintf(file,'| =========                 |                                                 |\n');
        fprintf(file,'| %s      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',mystring);
        fprintf(file,'|  %s    /   O peration     | Version:  v2012                                 |\n',mystring);
        fprintf(file,'|   %s  /    A nd           | Website:  www.openfoam.com                      |\n',mystring);
        fprintf(file,'|    %s/     M anipulation  |                                                 |\n',mystring);
        fprintf(file,'%s*---------------------------------------------------------------------------*/\n',mystring);
        fprintf(file,'FoamFile\n');
        fprintf(file,'{\n');
        fprintf(file,'    version     2.0;\n');
        fprintf(file,'    format      ascii;\n');
        fprintf(file,'    class       volVectorField;\n');
        fprintf(file,'    location    "3";\n');
        fprintf(file,'    object      traction;\n');
        fprintf(file,'}\n');
        fprintf(file,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
        fprintf(file,'dimensions      [1 -1 -2 0 0 0 0]; \n');
        fprintf(file,'internalField   uniform (0 0 0);\n \n');

        fprintf(file,'boundaryField\n');
        fprintf(file,'{\n');
        fprintf(file,'    free_Inside\n');
        fprintf(file,'    {\n');
        fprintf(file,'        type            calculated;\n');
        fprintf(file,'        value           uniform (0 0 0);\n');
        fprintf(file,'    }\n');
        fprintf(file,'    Thermowell\n');
        fprintf(file,'    {\n');
        fprintf(file,'        type            calculated;\n');
        fprintf(file,'        value           nonuniform List<vector> \n');

        fprintf(file,num2str(length(Traction)));
        fprintf(file,'\n( \n');

        for i = 1:length(Traction)
            fprintf(file, append('(',num2str(Traction(i,:)) ,')','\n' ));
        end

        fprintf(file, ') \n');
        fprintf(file, ';\n');
        fprintf(file, '    }\n');
        fprintf(file, '    fixed_Bottom \n');
        fprintf(file, '    { \n');
        fprintf(file, '        type         calculated\n');
        fprintf(file, '        value           uniform (0 0 0);\n');
        fprintf(file, '    } \n} \n');
        fprintf(file, '// ************************************************************************* //');

        fclose(file);
    end

end
