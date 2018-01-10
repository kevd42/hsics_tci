function [ cubename, lines,samples,bands,defaultBands, interleave,wavelengths] = parsecube( cubename )
    %allow either cube .raw or .hdr as input filename
    k = findstr('.hdr', cubename);
    if (length(k)>0) % hdr as input
        hdrname=cubename;
        cubename=regexprep(cubename,'\.hdr', '.raw');
    else
        cubename = cubename;
        hdrname=regexprep(cubename,'\.raw', '.hdr');
    end
    defaultBands=[49, 100, 150];

    fid = fopen(hdrname);
    tline = fgetl(fid);
    while ischar(tline)
        %disp(tline)
           M = sscanf(tline, 'samples = %d');
           if (length(M)>0)
               samples = M;
           end
           M = sscanf(tline, 'lines = %d');
           if (length(M)>0)
               lines = M;
           end
           M = sscanf(tline, 'bands = %d');
           if (length(M)>0)
               bands = M;
           end
           M = sscanf(tline, 'default bands = { %d, %d, %d} ');
           if (length(M)==3)
               defaultBands = M;
           end
           M = sscanf(tline, 'default bands = { %d, %d , %d} ');
           if (length(M)>0)
               defaultBands = M;
           end  
           M = sscanf(tline, 'interleave = %s');
           if (length(M)>0)
               interleave = M;
           end  
           M = sscanf(tline, 'Wavelength = %c');
           if (length(M)>0)
              
               for b=1:bands 
                   tline = fgetl(fid);
                   wl = sscanf(tline, '%f,');
                   wavelengths(b)=wl;
               end
           end  
        tline = fgetl(fid);
    end
    fclose(fid);
    %disp([ num2str(lines) 'x' num2str(samples) 'x' num2str(bands) ' cube'])

end

