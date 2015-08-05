classdef decompressBamFile < handle
    %DECOMPRESSBAMFILE object for decompressing BGZF blocks in a BAM file
    %
    %   OBJ = DECOMPRESSBAMFILE(filename) initializes a DECOMPRESSBAMFILE
    %   object for reading a BAM file.  OBJ stores the file identifiers to the
    %   BAM file and decompressed data blocks.  The DECOMPRESSBAMFILE class
    %   contains methods to decompress a BGZF block, and clean up file handles
    %   and temporary files upon deletion.
    %
    %   See also BAMREAD, BAMINFO
    
    %   Copyright 2009-2010 The MathWorks, Inc.

    
    properties
        bam_fid = -1 %File identifier for bam file
        BGZF_offset_end = inf %End read at this position in bam file
        data_offset_end %End read at this position in uncompressed data block
        is_eof %Reached the end of the bam file?
        
        tempfilename1 %File that stores BGZF block for decompression by gunzip
        tempfilename2 %File that stores BGZF block for decompression by gunzip
        file1_fid = -1 %File identifier for data file (file storing BGZF block)
        file2_fid = -1 %File identifier for data file (file storing BGZF block)
        data_file_size %Number of bytes in uncompressed block file
        switchFiles = 0 %Flag for use in bamread
        
        %Options used by baminfo and bamread
        parse_header = 1 %Parse header text using saminfo
        parse_dictionary = 1 %Parse reference dictionary
        tofile = 0 %Write data to SAM file
        sam_fid %File identifier of SAM file
    end
    
    methods
        
        function obj = decompressBamFile(filename)
            obj.bam_fid = fopen(filename, 'r', 'l');
            fread(obj.bam_fid, 1, '*uint8'); %Prepare file (because we usually read one extra byte to check for EOF)

            obj.tempfilename1 = tempname;
            obj.tempfilename2 = tempname;
        end
        
        function [fid last_block BGZF_file_pos_last] = decompressBlock(obj, bytesNeeded)
            %Previous position in BGZF block file, may be used for creating
            %index file
            BGZF_file_pos_last = ftell(obj.bam_fid)-1;
            try
                status = fseek(obj.bam_fid, 15, 'cof');
                if status
                    error('bioinfo:decompressBamFile:FseekFailed', ...
                        'FSEEK failed.  File may be incomplete.')
                end
                
                BSIZE = fread(obj.bam_fid, 1, 'uint16=>double');
                status = fseek(obj.bam_fid, -18, 'cof');
                if status
                    error('bioinfo:decompressBamFile:FseekFailed', ...
                        'FSEEK failed.')
                end
                data = fread(obj.bam_fid, BSIZE+1, '*uint8');
                
                one_more = fread(obj.bam_fid, 1, '*uint8'); %#ok<NASGU>
                last_block = feof(obj.bam_fid);
                obj.is_eof = feof(obj.bam_fid);
                
                fid2 = fopen([obj.tempfilename2 '.bam'], 'w', 'l');
                fwrite(fid2, data);
                fclose(fid2);
                
                gunzip([obj.tempfilename2 '.bam']);
                tmp = dir(obj.tempfilename2);
                obj.data_file_size = tmp.bytes;
                
                if obj.file2_fid ~= obj.file1_fid
                    fclose(obj.file2_fid);
                end
                obj.file2_fid = fopen(obj.tempfilename2, 'r', 'l');
                
                %Check whether we have reached BGZF_offset_end
                if ftell(obj.bam_fid)-1 > obj.BGZF_offset_end
                    last_block = 1;
                elseif ftell(obj.bam_fid) == obj.BGZF_offset_end
                    %We need this condition because when BGZF_offset_end is the
                    %last byte in the BAM file, data_offset_end = 0 and points to
                    %the start of a nonexistent decompressed block, instead of
                    %the end of the last block in the file.
                    obj.data_offset_end = obj.data_file_size;
                end
                
                if nargin > 1
                    if bytesNeeded < 0
                        error(message('bioinfo:decompressBamFile:NegativeBytesNeeded'))
                    end
                    copyData(obj, bytesNeeded);
                    fid = obj.file1_fid;
                else
                    fid = doSwitchFiles(obj);
                end
                
            catch ME
                error(message('bioinfo:decompressBamFile:DecompressionFailed'))
            end
        end
        
        function fid = doSwitchFiles(obj)
            if obj.file1_fid > 0
                obj.file1_fid = fclose(obj.file1_fid);
            end
            [obj.tempfilename1 obj.tempfilename2] = deal(obj.tempfilename2, obj.tempfilename1);
            obj.file1_fid = obj.file2_fid;
            fid = obj.file1_fid;
            obj.switchFiles = 0;
        end
        
        function moveFilePointer(obj, pos)
            status = fseek(obj.bam_fid, pos, 'bof');
            if status
                error('bioinfo:decompressBamFile:FseekFailed', ...
                    ['Failed to move pointer to position %d in file. ' ferror(obj.bam_fid)], pos)
            end
        end
        
        function delete(obj)
            if obj.bam_fid > 0
                fclose(obj.bam_fid);
            end
            if obj.file1_fid > 0
                fclose(obj.file1_fid);
            end
            if obj.file2_fid > 0 && obj.file2_fid ~= obj.file1_fid
                fclose(obj.file2_fid);
            end
            if exist(obj.tempfilename1, 'file')
                delete(obj.tempfilename1);
            end
            if exist([obj.tempfilename1 '.bam'], 'file')
                delete([obj.tempfilename1 '.bam']);
            end
            if exist(obj.tempfilename2, 'file')
                delete(obj.tempfilename2);
            end
            if exist([obj.tempfilename2 '.bam'], 'file')
                delete([obj.tempfilename2 '.bam']);
            end
            if ~isempty(obj.sam_fid)
                fclose(obj.sam_fid);
            end
        end
        
    end
    
    methods (Access = private)
        function copyData(obj, bytesNeeded)
            fid = fopen(obj.tempfilename1, 'a', 'l');
            %Read in remaining block and append to file 1
            tmp_data = fread(obj.file2_fid, bytesNeeded, '*uint8');
            fwrite(fid, tmp_data, '*uint8');
            fclose(fid);
            obj.switchFiles = 1; %Switch files after parsing this block
        end
    end
end
