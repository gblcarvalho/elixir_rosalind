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
    Probability.mendels_first_law(k, m, n) |> Float.to_string([compact: true])
  end

  def fib(args) do
    [n, k] = for i <- String.split(args), do: String.to_integer(i)
    DynamicProgramming.rabbits_recurrence_relations(n, k) |> Integer.to_string
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
  end
  
  def gc_content(dna) do
    counting_chars(dna, ["G","C"])
    |> Map.values
    |> Enum.reduce(0, &(&1 + &2))
    |> (fn(gc) -> (gc * 100)/ String.length(dna) end).()
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
    dnas = %{}
    content = String.split(content, ">", trim: true)

    Enum.reduce(content, dnas, fn(text, acc) ->
      [label | dna] = String.split(text, ["\r\n","\n"])
      dna = Enum.join(dna, "")
      Map.put(acc, label, dna)
    end)
  end
end


case File.read("test.txt") do
  {:ok, body}      -> IO.puts inspect (body |> String.strip |> FASTA.group)
  {:error, reason} -> IO.puts "Error when open the file: #{reason}"
end


#{:ok, file} = File.open "result.txt", [:write]

#case File.read("rosalind_fib.txt") do
##case File.read("test.txt") do
#  {:ok, body}      -> IO.binwrite file, (body |> String.strip |> Rosalind.fib)
#  {:error, reason} -> IO.puts "Error when open the file: #{reason}"
#end

#File.close file