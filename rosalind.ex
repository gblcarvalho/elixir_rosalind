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


{:ok, file} = File.open "result.txt", [:write]

case File.read("rosalind_revc.txt") do
#case File.read("test.txt") do
  {:ok, body}      -> IO.binwrite file, (body |> String.strip |> Rosalind.revc)
  {:error, reason} -> IO.puts "Error when open the file: #{reason}"
end

File.close file