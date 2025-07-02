#m2matrix.jl: parses CIGAR and MD tags for each alignment and populates correlated mutation matrix
using XAM
using StatsBase
using DataFrames
using ColorSchemes
using Gadfly
using ArgParse
using DelimitedFiles
using LinearAlgebra
using Base.Threads

#Using fixed-length buffers to avoid heap allocation
const MAX_READ_LENGTH = 10000
const MAX_BUFFER = 1000
const BAM_CIGAR_OPCODE = "MIDNSHP=X"
const BAM_MD_OPCODE = "=ACMGRSVTWYHKDBN"

#Structs to hold preallocated buffers
mutable struct CigarOp
    op::Char
    length::Int

    function CigarOp()
        op = Char(0)
        length = 0

        new(op, length)
    end
end

mutable struct CigarBuffer
    cigar_data::Vector{CigarOp}
    cigar_length::Int

    function CigarBuffer()
        cigar_data = Vector{CigarOp}(undef, MAX_BUFFER)
        for i in 1:MAX_BUFFER
            cigar_data[i] = CigarOp()
        end

        cigar_length = 0

        new(cigar_data, cigar_length)
    end
end

mutable struct MDOp
    op::Char
    length::Int
    seq::Vector{Char}
    seqlen::Int

    function MDOp()
        op = Char(0)
        length = 0
        seq = fill(Char(0), MAX_BUFFER)
        seqlen = 0

        new(op, length, seq, seqlen)
    end
end

mutable struct MDBuffer
    md_data::Vector{MDOp}
    md_length::Int
    current_field::Vector{Char}
    md_string::Vector{Char}
    md_string_length::Int

    function MDBuffer()
        md_data = Vector{MDOp}(undef, MAX_BUFFER)
        for i in 1:MAX_BUFFER
            md_data[i] = MDOp()
        end

        md_length = 0
        current_field = fill(Char(0), MAX_BUFFER)
        md_string = fill(Char(0), MAX_READ_LENGTH)
        md_string_length = 0

        new(md_data, md_length, current_field, md_string, md_string_length)
    end
end

mutable struct Mutation
    left::Int
    right::Int
    length::Int
    seq::Vector{Char}
    qual::Vector{Char}
    muttype::Char

    function Mutation()
        left = 0
        right = 0
        length = 0
        seq = fill(Char(0), MAX_BUFFER)
        qual = fill(Char(0), MAX_BUFFER)
        muttype = Char(0)

        new(left, right, length, seq, qual, muttype)
    end
end

mutable struct MutationBuffer
    mutation_data::Vector{Mutation}
    mutation_length::Int
    tmp_target_seq::Vector{Char}
    tmp_target_qual::Vector{Char}
    tmp_query_seq::Vector{Char}
    tmp_query_qual::Vector{Char}
    target_seq_overlap::Vector{Char}
    target_qual_overlap::Vector{Char}
    query_seq_overlap::Vector{Char}
    query_qual_overlap::Vector{Char}
    seq::Vector{Char}
    qual::Vector{Char}
    seqlen::Int
    mutation_pos::Vector{Int}
    mutation_pos_count::Int
    string_repr::Vector{Char}
    string_tmp::Vector{Char}
    string_repr_len::Int

    function MutationBuffer()
        mutation_data = Vector{Mutation}(undef, MAX_BUFFER)
        for i in 1:MAX_BUFFER
            mutation_data[i] = Mutation()
        end

        mutation_length = 0
        tmp_target_seq = fill(Char(0), MAX_BUFFER)
        tmp_target_qual = fill(Char(0), MAX_BUFFER)
        tmp_query_seq = fill(Char(0), MAX_BUFFER)
        tmp_query_qual = fill(Char(0), MAX_BUFFER)
        target_seq_overlap = fill(Char(0), MAX_BUFFER)
        target_qual_overlap = fill(Char(0), MAX_BUFFER)
        query_seq_overlap = fill(Char(0), MAX_BUFFER)
        query_qual_overlap = fill(Char(0), MAX_BUFFER)
        seq = fill(Char(0), MAX_READ_LENGTH)
        qual = fill(Char(0), MAX_READ_LENGTH)
        seqlen = 0
        mutation_pos = zeros(Int, MAX_BUFFER)
        mutation_pos_count = 0
        string_repr = fill(Char(0), MAX_BUFFER)
        string_tmp = fill(Char(0), MAX_BUFFER)
        string_repr_len = 0

        new(mutation_data, mutation_length, tmp_target_seq, tmp_target_qual, tmp_query_seq, tmp_query_qual,
            target_seq_overlap, target_qual_overlap, query_seq_overlap, query_qual_overlap, seq, qual, 
            seqlen, mutation_pos, mutation_pos_count, string_repr, string_tmp, string_repr_len)
    end
end

mutable struct MutationStats
    nm_dist::Matrix{Int}
    read_num::Vector{Int}
    muttype_num::Matrix{Int}
    read_num_pf::Vector{Int}
    read_num_pf_dist::Matrix{Int}

    function MutationStats(n_references::Int, max_reflen::Int)
        nm_dist = zeros(Int, n_references, MAX_READ_LENGTH)
        read_num = zeros(Int, n_references)
        muttype_num = zeros(Int, n_references, 3)
        read_num_pf = zeros(Int, n_references)
        read_num_pf_dist = zeros(Int, n_references, MAX_READ_LENGTH)

        new(nm_dist, read_num, muttype_num, read_num_pf, read_num_pf_dist)
    end
end

function SumMutationStats!(final_stats::MutationStats, stats::MutationStats)
    final_stats.nm_dist .+= stats.nm_dist
    final_stats.read_num .+= stats.read_num
    final_stats.muttype_num .+= stats.muttype_num
    final_stats.read_num_pf .+= stats.read_num_pf
    final_stats.read_num_pf_dist .+= stats.read_num_pf_dist
end

mutable struct ResultBuffer
    mutation_matrix_all::Array{Int}
    mutation_vector_all::Array{Int}
    effective_coverage_all::Matrix{Int}
    reflens::Vector{Int}
    refnames::Vector{String}
    n_references::Int
    max_reflen::Int
    mutation_stats::MutationStats

    function ResultBuffer(n_references::Int, max_reflen::Int, refnames::Vector{String}, reflens::Vector{Int})
        mutation_matrix_all = zeros(Int, n_references, max_reflen, max_reflen)
        mutation_vector_all = zeros(Int, n_references, 4, max_reflen)
        mutation_stats = MutationStats(n_references, max_reflen)
        effective_coverage_all = zeros(Int, n_references, max_reflen)
        
        new(mutation_matrix_all, mutation_vector_all, effective_coverage_all, reflens, refnames, n_references, max_reflen, mutation_stats)
    end
end

function SumResultBuffer!(final_result::ResultBuffer, result::ResultBuffer)
    final_result.mutation_matrix_all .+= result.mutation_matrix_all
    final_result.mutation_vector_all .+= result.mutation_vector_all
    final_result.effective_coverage_all .+= result.effective_coverage_all
    SumMutationStats!(final_result.mutation_stats, result.mutation_stats)
end

#Parse MD tag (outputs in-place to md_buffer)
function ParseMD!(md_buffer::MDBuffer)
    #First get references to preallocated md_buffer object
    md_data = md_buffer.md_data
    md_length = 0
    current_field = md_buffer.current_field
    md_string = md_buffer.md_string
    md_string_length = md_buffer.md_string_length #This is read-only

    last_char_type = 0 #Initially 0; 1 = number; 2 = not number
    current_field_length = 0
        
    #Breaks up MD string to vector of seq length and sequence
    for i in 1:md_string_length
        c = md_string[i]
        char_type = isdigit(c) ? 1 : 2
        if char_type == last_char_type #If current char type is last char type
            current_field_length += 1
            current_field[current_field_length] = c #Append character to current field
        else #If it changed, add to MDOp vector and reset current_field
            if current_field_length > 0 #If current field is not empty
                if isdigit(current_field[1]) #If starting with number, it's a match
                    seqlen = 0
                    for j in 1:current_field_length
                        digit = current_field_length - j
                        seqlen += parse(Int, current_field[j])*(10^digit) #change character to number
                    end
                    if seqlen != 0
                        md_length += 1
                        md_data[md_length].op = 'M'
                        md_data[md_length].length = seqlen
                        md_data[md_length].seqlen = 0
                    end
                elseif current_field[1] == '^' #If starting with ^, it's deletion
                    md_length += 1
                    md_data[md_length].op = 'D'
                    md_data[md_length].length = current_field_length - 1
                    md_data[md_length].seqlen = current_field_length - 1
                    for k in 2:current_field_length
                        md_data[md_length].seq[k-1] = current_field[k]
                    end
                else #Otherwise, it should be starting with A/T/C/G and it's substitution/insertion
                    md_length += 1
                    md_data[md_length].op = 'S'
                    md_data[md_length].length = current_field_length
                    md_data[md_length].seqlen = current_field_length
                    for k in 1:current_field_length
                        md_data[md_length].seq[k] = current_field[k]
                    end                 
                end
            end
            current_field_length = 1
            current_field[1] = c #Reset current field to current character
            last_char_type = char_type #Set last char type to current char type
        end
    end
    
    #Add the last MDOp entry from current_field
    if current_field_length > 0 #If current field is not empty
        if isdigit(current_field[1]) #If starting with number, it's a match
            seqlen = 0
            for j in 1:current_field_length
                digit = current_field_length - j
                seqlen += parse(Int, current_field[j])*(10^digit) #change character to number
            end
            if seqlen != 0
                md_length += 1
                md_data[md_length].op = 'M'
                md_data[md_length].length = seqlen
                md_data[md_length].seqlen = 0
            end
        elseif current_field[1] == '^' #If starting with ^, it's deletion
            md_length += 1
            md_data[md_length].op = 'D'
            md_data[md_length].length = current_field_length - 1
            md_data[md_length].seqlen = current_field_length - 1
            for k in 2:current_field_length
                md_data[md_length].seq[k-1] = current_field[k]
            end
        else #Otherwise, it should be starting with A/T/C/G and it's substitution/insertion
            md_length += 1
            md_data[md_length].op = 'S'
            md_data[md_length].length = current_field_length
            md_data[md_length].seqlen = current_field_length
            for k in 1:current_field_length
                md_data[md_length].seq[k] = current_field[k]
            end                 
        end
    end
    
    md_buffer.md_length = md_length
end

#The main parser function: it returns reference-centric mutation information from CIGAR and MD string
#Bug fixed and simplified implementation of parser in ShapeMapper:
#https://github.com/Weeks-UNC/shapemapper2/blob/master/internals/cpp-src/src/MutationParser.cpp
function ParseMutation!(mutation_buffer::MutationBuffer, cigar_buffer::CigarBuffer, md_buffer::MDBuffer, offset::Int, debug=0)
    #References to pre-allocated mutation_buffer, md_buffer, cigar_buffer
    mutation_data = mutation_buffer.mutation_data
    tmp_target_seq = mutation_buffer.tmp_target_seq
    tmp_target_qual = mutation_buffer.tmp_target_qual
    tmp_query_seq = mutation_buffer.tmp_query_seq
    tmp_query_qual = mutation_buffer.tmp_query_qual
    target_seq_overlap = mutation_buffer.target_seq_overlap
    target_qual_overlap = mutation_buffer.target_qual_overlap
    query_seq_overlap = mutation_buffer.query_seq_overlap
    query_qual_overlap = mutation_buffer.query_qual_overlap
    seq = mutation_buffer.seq
    qual = mutation_buffer.qual
    cigar_data = cigar_buffer.cigar_data
    cigar_length = cigar_buffer.cigar_length
    md_data = md_buffer.md_data
    md_length = md_buffer.md_length
    
    #Various indices and lengths to keep track
    tmp_op = -1
    tmp_length = 0
    tmp_mdop_index = -1
    tmp_seq_length = 0
    tmp_target_seq_length = 0
    mutation_length = 0
    in_match = false
    remaining_cigarop_length = 0
    query_index = 1
    mdop_index = 1
    target_index = 1
      
    #Iterate MD Ops until MD Op length adds up to CIGAR Op length. 
    #CIGAR alignment match, not sequence match, so need to keep track!
    for cigarop_index in 1:cigar_length
        if cigar_data[cigarop_index].op == 'M' #Process CIGAR Matches
            if !in_match #Reset when CIGAR and MD ops sum to same length
                in_match = true
                #Various stuff to hold for handling each MD within the same CIGAR op
                remaining_cigarop_length = 0
                tmp_op = md_data[mdop_index].op
                tmp_length = md_data[mdop_index].length

                if md_data[mdop_index].seqlen > 0
                    tmp_target_seq_length = md_data[mdop_index].seqlen
                    for i in 1:tmp_target_seq_length
                        tmp_target_seq[i] = md_data[mdop_index].seq[i]
                    end
                else
                    tmp_target_seq_length = 0
                end

                for i in 1:tmp_length
                    tmp_target_qual[i] = qual[query_index+i-1]
                    tmp_query_seq[i] = seq[query_index+i-1]
                    tmp_query_qual[i] = qual[query_index+i-1]
                end
                tmp_mdop_index = mdop_index
            end
            
            #Begin processing MD's for current CIGAR op
            remaining_cigarop_length += cigar_data[cigarop_index].length
            while (mdop_index <= md_length) && (remaining_cigarop_length > 0)
                if mdop_index != tmp_mdop_index #For each new MD op
                    tmp_op = md_data[mdop_index].op 
                    tmp_length = md_data[mdop_index].length

                    if md_data[mdop_index].seqlen > 0
                        tmp_target_seq_length = md_data[mdop_index].seqlen
                        for i in 1:tmp_target_seq_length
                            tmp_target_seq[i] = md_data[mdop_index].seq[i]
                        end
                    else
                        tmp_target_seq_length = 0
                    end
                    
                    for i in 1:tmp_length
                        tmp_target_qual[i] = qual[query_index+i-1]
                        tmp_query_seq[i] = seq[query_index+i-1]
                        tmp_query_qual[i] = qual[query_index+i-1]
                    end
                    tmp_mdop_index = mdop_index
                end

                overlap_length = 0
                if tmp_length > remaining_cigarop_length #MD op is longer than current CIGAR M op length
                    overlap_length = remaining_cigarop_length 
                else
                    overlap_length = tmp_length
                    mdop_index += 1
                end

                if tmp_op == 'M' #MD op is match
                    tmp_length = tmp_length - overlap_length
                else #Mismatch MD op
                    for i in 1:overlap_length
                        target_seq_overlap[i] = tmp_target_seq[i]
                        target_qual_overlap[i] = tmp_target_qual[i]
                        query_seq_overlap[i] = tmp_query_seq[i]
                        query_qual_overlap[i] = tmp_query_qual[i]
                    end
                    
                    for i in 1:(tmp_length-overlap_length)
                        tmp_target_seq[i] = tmp_target_seq[overlap_length+1]
                        tmp_target_qual[i] = tmp_target_qual[overlap_length+1]
                        tmp_query_seq[i] = tmp_query_seq[overlap_length+1]
                        tmp_query_qual[i] = tmp_query_qual[overlap_length+1]
                    end
                    
                    mutation_length += 1
                    mutation_data[mutation_length].left = offset+target_index
                    mutation_data[mutation_length].right = offset+target_index+overlap_length-1
                    mutation_data[mutation_length].length = overlap_length
                    for i in 1:overlap_length
                        mutation_data[mutation_length].seq[i] = query_seq_overlap[i]
                        mutation_data[mutation_length].qual[i] = query_qual_overlap[i]
                    end
                    mutation_data[mutation_length].muttype = 'S'
                end

                target_index += overlap_length
                query_index += overlap_length
                remaining_cigarop_length = remaining_cigarop_length - overlap_length
            end

            if (remaining_cigarop_length == 0) && (tmp_length == 0)
                in_match = false #CIGAR matches and MD tag sums to same length, flag to move onto next CIGAR
            end
        elseif cigar_data[cigarop_index].op == 'I' #Process CIGAR I's
            mutation_length += 1
            mutation_data[mutation_length].left = offset+target_index
            mutation_data[mutation_length].right = offset+target_index
            mutation_data[mutation_length].length = cigar_data[cigarop_index].length
            for i in 1:mutation_data[mutation_length].length
                mutation_data[mutation_length].seq[i] = seq[i+query_index-1]
                mutation_data[mutation_length].qual[i] = qual[i+query_index-1]
            end
            mutation_data[mutation_length].muttype = 'I'
            query_index += cigar_data[cigarop_index].length
        elseif cigar_data[cigarop_index].op == 'D' #Process CIGAR D's
            mutation_length += 1
            mutation_data[mutation_length].left = offset+target_index
            mutation_data[mutation_length].right = offset+target_index+cigar_data[cigarop_index].length-1
            mutation_data[mutation_length].length = 0
            mutation_data[mutation_length].muttype = 'D'
            target_index = target_index + cigar_data[cigarop_index].length
            mdop_index += 1
        elseif cigar_data[cigarop_index].op == 'N' #Process CIGAR N's
            target_index += 1             
        elseif cigar_data[cigarop_index].op == 'S' #Process CIGAR S's
            query_index += cigar_data[cigarop_index].length
        elseif cigar_data[cigarop_index].op == 'P' #Process CIGAR P's
            target_index += cigar_data[cigarop_index].length
        elseif cigar_data[cigarop_index].op == '=' #Process CIGAR ='s
            mdop_index += 1
        elseif cigar_data[cigarop_index].op == 'X' #Process CIGAR X's
            mutation_length += 1
            mutation_data[mutation_length].left = offset+target_index
            mutation_data[mutation_length].right = offset+target_index+cigar_data[cigarop_index].length-1
            mutation_data[mutation_length].length = cigar_data[cigarop_index].length - 1
            for i in 1:mutation_data[mutation_length].length
                mutation_data[mutation_length].seq[i] = md_data[mdop_index].seq[i]
                mutation_data[mutation_length].qual[i] = qual[i+query_index-1]
            end
            mutation_data[mutation_length].muttype = 'X'
            query_index += cigar_data[cigarop_index].length
            target_index += cigar_data[cigarop_index].length
            mdop_index += 1
        end
    end
    
    mutation_buffer.mutation_length = mutation_length
end

function ProcessRead!(
    result_buffer::ResultBuffer,
    mutation_buffer::MutationBuffer, 
    cigar_buffer::CigarBuffer, 
    md_buffer::MDBuffer, 
    record::BAM.Record, 
    min_mapq::Int, min_phred::Int, 
    min_mut_vector::Vector{Int}, 
    min_length_vector::Vector{Int},
    dump_file_handles::Vector{Vector{IOStream}},
    thread_id::Int) 
    #mutation_matrix_all::Array{Int})
    
    mutation_matrix_all = result_buffer.mutation_matrix_all
    mutation_vector_all = result_buffer.mutation_vector_all
    effective_coverage_all = result_buffer.effective_coverage_all
    mutation_stats = result_buffer.mutation_stats

    #mapq = UInt8((record.bin_mq_nl >> 8) & 0xff)
    mapq = record.mapq
    if (mapq >= min_mapq)
        #seqname_length = record.bin_mq_nl & 0xff #record.data starts at read name
        seqname_length = record.l_read_name #record.data starts at read name
        reference_id = record.refid + 1 #Ref ids are zero based
        start_pos = record.pos + 1 #zero-based leftmost coordinate
        seqlength = record.l_seq % Int #Long of read (not alignment)
        #n_cigar_ops = record.flag_nc & 0xFFFF #Total number of CIGAR ops
        n_cigar_ops = record.n_cigar_op #Total number of CIGAR ops

        #Various BAM data field offsets 
        #http://samtools.github.io/hts-specs/SAMv1.pdf section 4.2
        #record.data is byte array; these are index offsets
        data_size = record.block_size - 36 + sizeof(record.block_size) #total length of record.data; size of fixed-length fields add up to 36
        cigar_offset = seqname_length #following read name
        seq_offset = cigar_offset + n_cigar_ops * 4 #following cigar - uint32 per each cigar op
        qual_offset = seq_offset + cld(seqlength, 2) #following sequence - one byte for two bases (4 bit encoding per base)
        auxdata_offset = qual_offset + seqlength #following quality - one byte per base
        
        #Reinterpret CIGAR bytes and save to CigarOp struct vector in buffer
        align_length = 0
        cigar_buffer.cigar_length = n_cigar_ops
        for i in 1:n_cigar_ops
            cigar_bytes = reinterpret(UInt32, view(record.data, (cigar_offset+(i-1)*4+1):(cigar_offset+i*4)))[1] #4 bytes per each op; length=op_len<<4|op
            cigar_buffer.cigar_data[i].op = cigar_bytes % 4
            cigar_buffer.cigar_data[i].length = cigar_bytes >> 4
            
            #Does enumeration for mapping CIGAR op code; ‘MIDNSHP=X’→‘012345678’
            cigar_buffer.cigar_data[i].op = BAM_CIGAR_OPCODE[Int(cigar_buffer.cigar_data[i].op)+1]
            
            #Match and deletions add to alignment length
            if cigar_buffer.cigar_data[i].op == 'M'
                align_length += cigar_buffer.cigar_data[i].length
            elseif cigar_buffer.cigar_data[i].op == 'D'
                align_length += cigar_buffer.cigar_data[i].length
            end
        end
        
        if align_length >= min_length_vector[reference_id] #If alignment length (not read length) is longer than minimum       
            #Reinterpret sequence bytes and save to sequence Vector{Char} in mutation buffer
            for i in 1:(seqlength+1)
                if i % 2 == 1
                    mutation_buffer.seq[i] = record.data[i÷2+seq_offset] % 16
                else
                    mutation_buffer.seq[i] = record.data[i÷2+seq_offset] >> 4
                end
                if i >= 2
                    mutation_buffer.seq[i-1] = mutation_buffer.seq[i]    
                end
            end
            
            #Does enumeration for mapping base coding; ‘=ACMGRSVTWYHKDBN’→[0,15]        
            for i in 1:seqlength
                mutation_buffer.seq[i] = BAM_MD_OPCODE[Int(mutation_buffer.seq[i])+1]
            end

            #Reinterpret quality bytes and copy to quality Vector{Char} in mutation buffer
            for i in 1:seqlength
                mutation_buffer.qual[i] = record.data[i+qual_offset]
            end    
            
            mutation_buffer.seqlen = seqlength #Sequence length
            offset = Int(record.pos) #Start position

            #Get the MD tag from AUX data
            #AUX field format is documented in http://samtools.github.io/hts-specs/SAMv1.pdf section 4.2.4
            #Each field is two-character tag followed by one character specifying type, then value
            #Structure looks like this. tag (2 bytes) : type (1 byte) : data
            md_buffer.md_string_length = 0 #Clear MD string 
            md_offset = auxdata_offset + 1
            for i in 1:(data_size-auxdata_offset)
                if md_offset <= (data_size - 3)
                    found_md = false
                    found_nm = false
                    #We are looking for: data[pos]='M', data[pos+1]='D'            
                    if (record.data[md_offset] == 77) && (record.data[md_offset+1] == 68) 
                        found_md = true
                    end

                    if (record.data[md_offset] == 78) && (record.data[md_offset+1] == 77)
                        found_nm = true 
                    end

                    auxtype = Char(record.data[md_offset+2])
                    md_offset += 3 #Now we are at value
                    
                    #Aux data types; defined sizes
                    if auxtype == 'A'
                        md_offset += 1
                    end

                    if auxtype == 'c' || auxtype == 'C'
                        if found_nm
                            nm = reinterpret(UInt8, view(record.data, md_offset:md_offset))[1] #1 bytes
                            mutation_stats.nm_dist[reference_id, nm+1] += 1
                        end
                        md_offset += 1
                    end

                    if auxtype == 's' || auxtype == 'S'
                        if found_nm
                            nm = reinterpret(UInt16, view(record.data, md_offset:(md_offset+1)))[1] #2 bytes
                            mutation_stats.nm_dist[reference_id, nm+1] += 1
                        end
                        md_offset += 2
                    end

                    if auxtype == 'i' || auxtype == 'I'
                        if found_nm
                            nm = reinterpret(UInt32, view(record.data, md_offset:(md_offset+3)))[1] #4 bytes
                            mutation_stats.nm_dist[reference_id, nm+1] += 1
                        end
                        md_offset += 4
                    end

                    if auxtype == 'f'
                        md_offset += 4
                    end

                    if auxtype == 'd'
                        md_offset += 8
                    end
                    
                    #Now aux data types that are variable sizes                             
                    if auxtype == 'Z' || auxtype == 'H' #NULL-terminated string encoding
                        while record.data[md_offset] != 0x00  # NULL-terminated
                            if found_md                        
                                md_buffer.md_string_length += 1
                                md_buffer.md_string[md_buffer.md_string_length] = Char(record.data[md_offset])
                            end
                            md_offset += 1
                        end
                        md_offset += 1 #Termination byte. #TODO: double check this is correct when there are situations where multiple AUX arrays are presents
                    end
                    
                    if auxtype == 'B' #This means array
                        #tag (2 bytes) : B (1 byte) : element type (1 byte) : count (4 bytes) : array data
                        element_type = Char(record.data[md_offset]) #current offset is element type
                        element_size = element_type == 'c' || element_type == 'C' ? 1 : #uint8/int8
                                 element_type == 's' || element_type == 'S' ? 2 : #uint16/int16
                                 element_type == 'i' || element_type == 'I' || eltyp == 'f' ? 4 :#uint32/int32/float32
                                 error("invalid type tag: '$(Char(eltyp))'")
                        md_offset += 1 #Now we are at count
                        element_count = reinterpret(UInt32, view(record.data, md_offset:(md_offset+3)))[1] #4 bytes
                        md_offset += 4 #Now we are at data
                        md_offset += element_size * element_count #Now we are done
                    end

                    if found_md && found_nm
                        break #Don't need to keep going after MD is saved
                    end
                end
            end
            
            mutation_stats.read_num[reference_id] += 1

            ParseMD!(md_buffer)
            ParseMutation!(mutation_buffer, cigar_buffer, md_buffer, offset) #TODO: add stuff to catch errors before this

            ##Effective coverage
            view(effective_coverage_all, reference_id, start_pos:(start_pos+align_length-1)) .+= 1            
            
            mutation_data = mutation_buffer.mutation_data

            #Build a 1D vector of mutation positions
            mutation_buffer.mutation_pos_count = 0 #Clear mutation position vecetor
            for i in 1:mutation_buffer.mutation_length
                if mutation_data[i].muttype == 'S' #Substitution only
                    if mutation_data[i].seq[1] != 'N' #Ignore N's
                        if Int(mutation_data[i].qual[1]) >= min_phred #PHRED cutoff
                            mutation_buffer.mutation_pos_count += 1
                            mutation_buffer.mutation_pos[mutation_buffer.mutation_pos_count] = mutation_data[i].left

                            #Mutation 1D vector
                            if mutation_data[i].seq[1] == 'A'
                                mutation_vector_all[reference_id, 1, mutation_data[i].left] += 1
                            elseif mutation_data[i].seq[1] == 'T'
                                mutation_vector_all[reference_id, 2, mutation_data[i].left] += 1
                            elseif mutation_data[i].seq[1] == 'G'
                                mutation_vector_all[reference_id, 3, mutation_data[i].left] += 1
                            elseif mutation_data[i].seq[1] == 'C'
                                mutation_vector_all[reference_id, 4, mutation_data[i].left] += 1
                            end

                            mutation_stats.muttype_num[reference_id, 1] += 1
                        end
                    end
                elseif mutation_data[i].muttype == 'D'
                    mutation_stats.muttype_num[reference_id, 2] += 1
                elseif mutation_data[i].muttype == 'I'
                    mutation_stats.muttype_num[reference_id, 3] += 1
                end
            end
            
            #Update the 2D matrix by incrementing each cell in M2 matrix for all possible pairs of mutations
            mutation_stats.read_num_pf_dist[reference_id, mutation_buffer.mutation_pos_count+1] += 1
            if (mutation_buffer.mutation_pos_count >= min_mut_vector[reference_id]) && (mutation_buffer.mutation_pos_count >= 2)
                mutation_stats.read_num_pf[reference_id] += 1
                for i in 1:mutation_buffer.mutation_pos_count
                    for j in 1:mutation_buffer.mutation_pos_count
                        mutation_matrix_all[reference_id, mutation_buffer.mutation_pos[i], mutation_buffer.mutation_pos[j]] += 1
                    end
                end
            end

            #IR
            string_repr = mutation_buffer.string_repr
            string_tmp = mutation_buffer.string_tmp
            fill!(string_tmp, '.')

            #Sub and dels
            for i in 1:mutation_buffer.mutation_length
                if mutation_data[i].muttype == 'S' #Substitution only
                    if mutation_data[i].seq[1] != 'N' #Ignore N's
                        if Int(mutation_data[i].qual[1]) >= min_phred
                            string_tmp[mutation_data[i].left] = mutation_data[i].seq[1]   
                        end
                    end
                elseif mutation_data[i].muttype == 'D'
                    view(string_tmp, mutation_data[i].left:mutation_data[i].right) .= '-'
                end        
            end

            #Insertions
            ins_count = 0
            repr_pos = 1
            tmp_pos = 1
            repr_len = result_buffer.reflens[reference_id]
            ins_start_pos = 0
            ins_end_pos = 0
            for i in 1:mutation_buffer.mutation_length
                if mutation_data[i].muttype == 'I' #Substitution only
                    ins_start_pos = mutation_data[i].left + ins_count
                    ins_end_pos = ins_start_pos + (2 + mutation_data[i].length - 1)
                    view(string_repr, (ins_start_pos+1):(ins_end_pos-1)) .=  mutation_data[i].seq[1:mutation_data[i].length]
                    string_repr[ins_start_pos] = '['
                    string_repr[ins_end_pos] = ']'
                    view(string_repr, repr_pos:(ins_start_pos-1)) .= view(string_tmp, tmp_pos:(mutation_data[i].left-1))
                    tmp_pos = mutation_data[i].left
                    repr_pos = ins_end_pos + 1
                    ins_count += (2 + mutation_data[i].length)
                end
            end
            if ins_count == 0
                view(string_repr, 1:repr_len) .= view(string_tmp, 1:repr_len)
                mutation_buffer.string_repr_len = repr_len
            else
                mutation_buffer.string_repr_len = ins_count + repr_len
                view(string_repr, (ins_end_pos+1):(ins_count+repr_len)) .= view(string_tmp, (tmp_pos:repr_len))
            end   
            
            f = dump_file_handles[reference_id][thread_id]
            accepted = (mutation_buffer.mutation_pos_count >= min_mut_vector[reference_id]) && (mutation_buffer.mutation_pos_count >= 2)
            write(f, join(string_repr[1:mutation_buffer.string_repr_len]))
            write(f, ',')
            write(f,  string(accepted))
            write(f, '\n')
        end
    end
end

function rownorm(count_matrix::SubArray{Int})
    matrix_size = size(count_matrix)[1]
    diag_val = diag(count_matrix)
    if sum(count_matrix) == 0
        return zeros(Float64, matrix_size, matrix_size)#, sizefactors, lib_sizes
    end
    
    scale_matrix = reshape(repeat(diag_val, matrix_size), matrix_size, matrix_size) #Fill a tmp matrix
    mutation_matrix_norm = log2.(transpose((count_matrix .+ 0.5) ./ (scale_matrix .+ 0.5) .* 1e6)) #Thus the final scaling is: log2 (counts+0.5) / (scale) * 1e6

    return mutation_matrix_norm
end

function ParseCmdArgs()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--bam"
            help = "Input BAM file (uses CIGAR and MD tags to parse mutations); assumes merged or single-end reads"
            required = true
        "--out"
            help = "Output path prefix (make sure path exists)"
            required = true
        #"--ref"
        #    help = "List of reference names"
        #    required = true
        "--ltrim"
           help = "Trim left"
           arg_type = Int
           default = 0
        "--rtrim"
           help = "Trim right"
           arg_type = Int
           default = 0
        "--minmut"
            help = "Minumum number of mutations for a read to be used; expressed as fraction of reference length (this or at least 2)"
            arg_type = Float64
            default = 0.01
        "--minphred"
            help = "Minimum PHRED score for a base to be counted"
            arg_type = Int
            default = 30
        "--minmapq"
            help = "Minimum MAPQ score for a read to be used"
            arg_type = Int
            default = 0
        "--minlength"
            help = "Minimum alignment length for a read to be used; expressed as fraction of reference length"
            arg_type = Float64
            default = 0.5
        "--plot"
            help = "Output a plot of the heatmap"
            arg_type = Bool
            default = false
    end
    return parse_args(s)
end

function main()
    parsed_args = ParseCmdArgs()
    bam_file = parsed_args["bam"]
    min_mapq = parsed_args["minmapq"]
    min_phred = parsed_args["minphred"]
    min_mut = parsed_args["minmut"]
    min_length = parsed_args["minlength"]
    out_prefix = parsed_args["out"]
    plot_flag = parsed_args["plot"]
    #ref_file = parsed_args["ref"]
    ltrim = parsed_args["ltrim"]
    rtrim = parsed_args["rtrim"]
    
    #ref_out = readlines(ref_file)
    bam_reader = open(BAM.Reader, bam_file) #Reader from XAM
    
    #Make 3D array {references, reflen, reflen}. Careful with max reference size.
    n_references = length(bam_reader.refseqnames)
    max_reflen = maximum(bam_reader.refseqlens)
    mutation_matrix_final = zeros(Int, n_references, max_reflen, max_reflen)
    
    result_final = ResultBuffer(n_references, max_reflen, bam_reader.refseqnames, bam_reader.refseqlens)

    min_length_vector::Vector{Int} = [floor(min_length * len) for len in bam_reader.refseqlens] #Minimum lengths
    min_mut_vector::Vector{Int} = [floor(min_mut * len) for len in bam_reader.refseqlens] #Minimum mutations
    
    worker_threads = Threads.nthreads() - 1
    
    #TODO: deal with thread safe write later
    #IR dump file handles
    file_handles_threads = Vector{Vector{IOStream}}()
    for i in 1:n_references
        file_handles = Vector{IOStream}()
        for j in 1:worker_threads
            f = open(string(out_prefix, bam_reader.refseqnames[i], "_", string(j), ".dump"), "w")
            push!(file_handles, f)
        end
        push!(file_handles_threads, file_handles)
    end

    #Merging multithreading from April
    finished = false
    ready = [0 for _ = 1:worker_threads]
    buffer_size = 128 #TODO: take a closer look at this
    slots = [[BAM.Record() for _ = 1:buffer_size] for _ = 1:worker_threads]
    reader = Threads.@spawn begin
        while !eof(bam_reader)
            for i in 1:worker_threads
                if ready[i] == 0
                    count = 0
                    while count < buffer_size && !eof(bam_reader)
                        #Expect this to be bottlenecking at some point. BGZF decompression is happening here.
                        #Also this reader allocates a lot. Though total size isn't too bad, it does scale with read number.
                        #TODO: test multiple decompression threads. Also should write my own BAM reader that avoids allocation.
                        read!(bam_reader, slots[i][count + 1])
                        count += 1
                    end
                    ready[i] = count
                end
            end
            yield()
        end
        finished = true
    end
       
    sync_lock = ReentrantLock()
    Threads.@threads for i = 1:worker_threads
        mutation_matrix_all = zeros(Int, n_references, max_reflen, max_reflen)
        
        #Initialize the buffers. Just zero constructing bunch of vectors.
        mutation_buffer = MutationBuffer()
        md_buffer = MDBuffer()
        cigar_buffer = CigarBuffer()
        result_buffer = ResultBuffer(n_references, max_reflen, bam_reader.refseqnames, bam_reader.refseqlens)
        
        #Now do the work
        while !finished || ready[i] > 0
            available = ready[i]
            if available == 0
                yield()
                continue
            end
            for s in 1:available
                ProcessRead!(result_buffer, mutation_buffer, cigar_buffer, md_buffer, slots[i][s], min_mapq, min_phred, min_mut_vector, min_length_vector, file_handles_threads, i)
            end
            ready[i] = 0
        end

        #Lock and add to shared final matrix
        lock(sync_lock)
        #mutation_matrix_final .+= mutation_matrix_all
        SumResultBuffer!(result_final, result_buffer)
        unlock(sync_lock)
    end
    wait(reader)
   
    #Output M2 matrix
    header = ["ref_name", "mean_edit_dist", "mean_sub", "used_reads", "used_reads_pf", "total_sub", "total_del", "total_ins"]
    println(join(header, '\t'))
    for i in 1:n_references
        reference_name = bam_reader.refseqnames[i]

        #if reference_name in ref_out
        effective_length = bam_reader.refseqlens[i]
        left_pos = ltrim + 1
        right_pos = effective_length - rtrim
        mutation_matrix = view(result_final.mutation_matrix_all, i, left_pos:right_pos, left_pos:right_pos)
        mut1d = view(result_final.mutation_vector_all, i, :, left_pos:right_pos)
        coverage = view(result_final.effective_coverage_all, i, left_pos:right_pos)
        read_num_pf_dist = view(result_final.mutation_stats.read_num_pf_dist, i, 1:(effective_length+1))
        mut1d_all = sum(mut1d, dims=1)[1, :]
        reac1d =  mut1d_all ./ coverage
        
        #Output sub histogram
        outfile_subs = string(out_prefix, reference_name, "_subs.txt")
        writedlm(outfile_subs, read_num_pf_dist)

        #Output 1D vector
        outfile_mut1d = string(out_prefix, reference_name, "_mut1d.txt")
        writedlm(outfile_mut1d, mut1d)

        #Effective coverage
        outfile_coverage = string(out_prefix, reference_name, "_cov.txt")
        writedlm(outfile_coverage, coverage)

        #1D reactivity scaled by effective coverage
        outfile_reac = string(out_prefix, reference_name, "_reac.txt")
        writedlm(outfile_reac, reac1d)

        #Print read number and other mutation stats
        #Calculate mean of edit distances
        edit_distances = view(result_final.mutation_stats.nm_dist, i, 1:effective_length)
        sum_edit = 0
        for j in 1:effective_length
            sum_edit += (j-1) * edit_distances[j]
        end
        avg_edit = sum_edit / sum(edit_distances)

        subs = view(result_final.mutation_stats.read_num_pf_dist, i, 1:effective_length)
        sum_sub = 0
        for j in 1:effective_length
            sum_sub += (j-1) * subs[j]
        end
        avg_sub = sum_sub / sum(subs)

        summary_tx = [reference_name, 
            string(avg_edit),
            string(avg_sub),
            string(result_final.mutation_stats.read_num[i]),
            string(result_final.mutation_stats.read_num_pf[i]),
            string(result_final.mutation_stats.muttype_num[i, 1]),
            string(result_final.mutation_stats.muttype_num[i, 2]),
            string(result_final.mutation_stats.muttype_num[i, 3])]
        println(join(summary_tx, '\t'))

        if plot_flag
            #Plot substitution hist
            yend = effective_length
            for i in 1:effective_length
                if read_num_pf_dist[yend] != 0
                    break
                else
                    yend -= 1
                end
            end
            df = DataFrame(
                n_sub=1:yend, 
                value=read_num_pf_dist[1:yend])
            p = plot(df, x="n_sub", y="value", Geom.bar,
                Guide.xlabel("Number of substitutions"), Guide.ylabel("Count"),
                Guide.title(reference_name))
            outfile_sub_plot = string(out_prefix, reference_name, "_sub_plots.svg")
            img = SVG(outfile_sub_plot, 6inch, 4inch)
            draw(img, p)

            if sum(reac1d) > (effective_length)     
                #Plot 1D
                df = DataFrame(
                    pos=1:effective_length,
                    value=reac1d[1:effective_length])
                p = plot(df, x="pos", y="value", Geom.bar,
                    Guide.xlabel("Position"), Guide.ylabel("Reactivity"),
                    Guide.title(reference_name),
                    Scale.x_continuous(minvalue=1, maxvalue=effective_length) )
                outfile_reac_plot = string(out_prefix, reference_name, "_reac.svg")
                img = SVG(outfile_reac_plot, 8inch, 3inch)
                draw(img, p)                
            end

            #mutation_matrix_norm, sizefactors, lib_sizes = TMMNorm(mutation_matrix) #TMM scaling wrt column
            mutation_matrix_norm = rownorm(mutation_matrix) #scaling wrt row
            z_transform = StatsBase.fit(ZScoreTransform, mutation_matrix_norm, dims=2) #Z-scaling wrt col
            mutation_matrix_norm_z = StatsBase.transform(z_transform, mutation_matrix_norm)
            
            #Write the matrices as delimited text files
            outfile_counts = string(out_prefix, reference_name, ".mat")
            writedlm(outfile_counts, mutation_matrix)
            
            outfile_norm = string(out_prefix, reference_name, "_norm.mat")
            writedlm(outfile_norm, mutation_matrix_norm)         
            
            outfile_norm_z = string(out_prefix, reference_name, "_norm_z.mat")
            writedlm(outfile_norm_z, mutation_matrix_norm_z)
            
            if sum(mutation_matrix) > (effective_length ^ 2) 
                #Plot z-scaled matrix as heatmap
                matrix_size = size(mutation_matrix_norm_z)[1]
                df = DataFrame(
                    row=repeat(1:matrix_size, outer=matrix_size), 
                    col=repeat(1:matrix_size, inner=matrix_size),
                    value=vec(mutation_matrix_norm_z))
                cpalette(p) = get(ColorSchemes.magma, p)
                p = plot(df, x="col", y="row", color="value", Geom.rectbin,
                    Coord.cartesian(yflip=true),    
                    Scale.color_continuous(minvalue=-3, maxvalue=3, colormap=cpalette),
                    Guide.xlabel(""), Guide.ylabel(""), Guide.colorkey(""),
                    Guide.title(reference_name))
                outfile_norm_z_plot = string(out_prefix, reference_name, "_norm_z_plot.svg")
                img = SVG(outfile_norm_z_plot, 6inch, 6inch)
                draw(img, p)
            end
        end
    end
end

main()
