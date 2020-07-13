module Assignment07

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

function composition(seq)
    sequence = normalizeDNA(seq) # make uppercase string, check invalid bases
    a = c = g = t = n = 0 # sets all 4 variables to `0`
    for base in sequence
        a = count(==('A'), sequence)
        c = count(==('C'), sequence)
        g = count(==('G'), sequence)
        t = count(==('T'), sequence)
        n = count(==('N'), sequence)
    end
    composition = Dict('A' => a, 'G' => g, 'T' => t, 'N' => n, 'C' => c)
    return composition
end

function gc_content(seq)
    comp = composition(seq)
    gs = comp['G']
    cs = comp['C']
    seqlength = length(seq)
    return (gs+cs) / (seqlength)
end

function complement(base::Char)
    complements = Dict('A' => 'T',
                       'T' => 'A',
                       'G' => 'C',
                       'C' => 'G',
                       'N' => 'N')
    
    !(base in keys(complements)) && error("Invalid base $base")
    return complements[base]
end

function complement(seq::AbstractString)
    seq = normalizeDNA(seq)
    seq = collect(seq)
    complement_seq = []
    for base in eachindex(seq)
        new_seq = complement(seq[base])
        push!(complement_seq, new_seq)
    end
    complement_seq = join(complement_seq)
    return complement_seq
end

function reverse_complement(seq)
    seq = normalizeDNA(seq)
    reverse_sequence = reverse(seq)
    sequence = collect(reverse_sequence)
    reverse_complement = map(complement, sequence)
    reverse_complement = join(reverse_complement)
end

function parse_fasta(path)
    header = String[]
    sequence = []
    seqonly = []
    for line in eachline(path)
        if startswith(line, '>')
            if length(header) > length(sequence)
                seqonlynew = join(seqonly)
                push!(sequence, seqonlynew)
                seqonly = []
            end
            headerline = split(line, '>')
            headerline = headerline[2]
            push!(header, headerline)
        elseif line == ""
            continue
        else
            #if line is not a header or an empty string, push it to seqonly
            #keep pushing these lines until you reach another header or empty string line
            #then join all of these lines together and push the sequence to sequence
            line = normalizeDNA(line)
            push!(seqonly, line) 
        end
    end
    seqonlynew = join(seqonly)
    push!(sequence, seqonlynew)
    return(header, sequence)
end

end # module Assignment07