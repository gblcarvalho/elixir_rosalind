defmodule Rosalind do
  def dna(dna) do
    result = DNA.counting_acgt(dna)
    "#{result["A"]} #{result["C"]} #{result["G"]} #{result["T"]}"
  end

  def rna(dna) do
    DNA.transcribe_dna_in_rna(dna)
  end

  def revc(dna) do
    DNA.reverse_complement(dna)
  end

  def iprb(args) do
    [k, m, n] = for i <- String.split(args), do: String.to_integer(i)
    Probability.mendels_first_law(k, m, n) |> Float.to_string
  end

  def fib(args) do
    [n, k] = for i <- String.split(args), do: String.to_integer(i)
    DynamicProgramming.rabbits_recurrence_relations(n, k) |> Integer.to_string
  end

  def gc(args) do
    {max_label, max_gc} = 
	  args
      |> FASTA.group
      |> DNA.highest_gc_content

    "#{max_label}\n#{Float.to_string(max_gc)}"
  end

  def prot(rna) do
    DNA.encoding_rna_into_aminoacid(rna)
  end
  
  def subs(args) do
    [dna, motif] = String.split(args, ["\r\n","\n"])

    (for i <- DNA.finding_locations_motif(dna, motif), into: "", do: "#{i} ")
    |> String.strip
  end
  
  def hamm(args) do
    [dna_one, dna_two] = String.split(args, ["\r\n","\n"])
    DNA.hamming_distance(dna_one, dna_two)
	|> Integer.to_string
  end
  
  def prtm(protein) do
    protein
	|> DNA.protein_mass
	|> Float.to_string([decimals: 5, compact: true])
  end
  
  def perm(number) do
    final = String.to_integer(number)
	
	permutations = 
	  (for i <- 1..final, do: i)
	  |> Combinatorics.permutation
	
	number_permutations = 
	  permutations
	  |> length
	  |> Integer.to_string
	
	permutations = (for p <- permutations, do: Enum.join(p, " "))
	number_permutations <> "\n" <> Enum.join(permutations, "\n")
  end
end


defmodule DNA do
  def counting_acgt(dna) do
    counting_chars(dna, ["A","C","G","T"])
  end

  def transcribe_dna_in_rna(dna) do
    String.replace(dna, "T", "U")
  end

  def reverse_complement(dna) do
    complements = %{"A" => "T", "T" => "A", "C" => "G", "G" => "C"}
    complement = for n <- String.codepoints(dna), into: "", do: complements[n]
    String.reverse(complement)
  end

  def highest_gc_content(dnas) do
    (for {label, dna} <- dnas, do: {label, gc_content(dna)})
    |> Enum.max_by(fn(gc) -> elem(gc, 1) end)
  end

  def gc_content(dna) do
    counting_chars(dna, ["G","C"])
    |> Map.values
    |> Enum.reduce(0, &(&1 + &2))
    |> (fn(gc) -> (gc * 100)/ String.length(dna) end).()
  end

  def encoding_rna_into_aminoacid(rna) do
    rna_codon_table = 
	  %{"UUU" => "F", "CUU" => "L", "AUU" => "I", "GUU" => "V",
        "UUC" => "F", "CUC" => "L", "AUC" => "I", "GUC" => "V",
        "UUA" => "L", "CUA" => "L", "AUA" => "I", "GUA" => "V",
        "UUG" => "L", "CUG" => "L", "AUG" => "M", "GUG" => "V",
        "UCU" => "S", "CCU" => "P", "ACU" => "T", "GCU" => "A",
        "UCC" => "S", "CCC" => "P", "ACC" => "T", "GCC" => "A",
        "UCA" => "S", "CCA" => "P", "ACA" => "T", "GCA" => "A",
        "UCG" => "S", "CCG" => "P", "ACG" => "T", "GCG" => "A",
        "UAU" => "Y", "CAU" => "H", "AAU" => "N", "GAU" => "D",
        "UAC" => "Y", "CAC" => "H", "AAC" => "N", "GAC" => "D",
        "UAA" => "Stop", "CAA" => "Q", "AAA" => "K", "GAA" => "E",
        "UAG" => "Stop", "CAG" => "Q", "AAG" => "K", "GAG" => "E",
        "UGU" => "C", "CGU" => "R", "AGU" => "S", "GGU" => "G",
        "UGC" => "C", "CGC" => "R", "AGC" => "S", "GGC" => "G",
        "UGA" => "Stop", "CGA" => "R", "AGA" => "R", "GGA" => "G",
        "UGG" => "W", "CGG" => "R", "AGG" => "R", "GGG" => "G"}

    do_encoding_rna_into_aminoacid(rna, rna_codon_table, "")
  end
  
  def finding_locations_motif(dna, motif) do
    Regex.scan(~r{(?=#{motif})}, dna, return: :index, capture: :first)
    |> Enum.map(fn (x) -> 
         [{index, _} | _tail] = x
         index + 1
       end)
  end
  
  def hamming_distance(dna_one, dna_two) do
    Enum.zip(String.codepoints(dna_one), String.codepoints(dna_two))
    |> Enum.reduce(0, fn(x, acc) ->
	     {n_one, n_two} = x
	     cond do
		   n_one == n_two -> acc
		   true -> acc + 1
		 end
	   end)
  end
  
  def protein_mass(protein) do
    monoisotopic_mass =
      %{"A" => 71.03711,
		"C" => 103.00919,
		"D" => 115.02694,
		"E" => 129.04259,
		"F" => 147.06841,
		"G" => 57.02146,
		"H" => 137.05891,
		"I" => 113.08406,
		"K" => 128.09496,
		"L" => 113.08406,
		"M" => 131.04049,
		"N" => 114.04293,
		"P" => 97.05276,
		"Q" => 128.05858,
		"R" => 156.10111,
		"S" => 87.03203,
		"T" => 101.04768,
		"V" => 99.06841,
		"W" => 186.07931,
		"Y" => 163.06333}

    protein
	|> String.codepoints
	|> Enum.reduce(0, fn(x, acc) ->
	     acc + monoisotopic_mass[x]
	   end)
  end

  defp do_encoding_rna_into_aminoacid("", _table, protein), do: protein
  defp do_encoding_rna_into_aminoacid(rna, table, protein) do
    case table[String.slice(rna, 0, 3)] do
      "Stop" -> protein
      nil -> do_encoding_rna_into_aminoacid(String.slice(rna, 1..-1),
                                            table,
                                            protein <> String.slice(rna, 1, 1))
      p -> do_encoding_rna_into_aminoacid(String.slice(rna, 3..-1), table, protein <> p)
    end
  end

  defp counting_chars(string, chars) do
    result = for char <- chars, into: %{}, do: {char, 0}
    Enum.reduce(String.codepoints(string), result, fn(letter, acc) ->
      cond do
        Map.has_key?(acc, letter) -> Map.update(acc, letter, 1, &(&1 + 1))
        true -> acc
      end
    end)
  end

end


defmodule Combinatorics do
  def permutation([]), do: [[]]  
  def permutation(list) when is_list(list) do
    for head <- list, tail <- permutation(list -- [head]), do: [head | tail]
  end
end


defmodule Probability do
  def mendels_first_law(k, m, n) do
    total = k + m + n
    kp = k/total
    mp = ((m/total) * (k/(total-1))) + ((m/total) * ((m-1)/(total-1)) * 0.75) + ((m/total) * (n/(total-1)) * 0.5)
    np = ((n/total) * (k/(total-1))) + ((n/total) * (m/(total-1)) * 0.5)
    kp + mp + np
  end
end


defmodule DynamicProgramming do
  def rabbits_recurrence_relations(n, k) do
    case n do
      1 -> 1
      2 -> 1
      _ -> rabbits_recurrence_relations(n-1, k) +
           (rabbits_recurrence_relations(n-2, k) * k)
    end
  end
end


defmodule FASTA do
  def group(content) do
    content = String.split(content, ">", trim: true)

    Enum.reduce(content, [], fn(text, acc) ->
      [label | dna] = String.split(text, ["\r\n","\n"])
      dna = Enum.join(dna, "")
      acc ++ [{label, dna}]
    end)
  end
end


{:ok, file} = File.open "result.txt", [:write]

case File.read("rosalind_perm.txt") do
#case File.read("test.txt") do
  {:ok, body}      -> IO.binwrite file, (body |> String.strip |> Rosalind.perm)
  {:error, reason} -> IO.puts "Error when open the file: #{reason}"
end

File.close file