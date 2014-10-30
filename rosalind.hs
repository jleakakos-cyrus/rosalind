import Data.List
import Data.List.Split

normalize :: [Int] -> String
normalize = intercalate " " . map show

-- [A, C, G, T]
-- Counts of nucleobases
dna :: String -> [Int]
dna = foldr (\letter acc -> update letter acc) [0,0,0,0]
  where update 'A' [a, c, g, t] = [a+1, c, g, t]
        update 'C' [a, c, g, t] = [a, c+1, g, t]
        update 'G' [a, c, g, t] = [a, c, g+1, t]
        update 'T' [a, c, g, t] = [a, c, g, t+1]
        update _   [a, c, g, t] = [a, c, g, t]
                                
-- T replaced with U
-- DNA to RNA
rna :: String -> String
rna = map (\letter -> if letter == 'T' then 'U' else letter)

-- reverse string, A <-> T, C <-> G
-- Compliment a strand of DNA
revc :: String -> String
revc = map update . reverse
  where update 'A' = 'T'
        update 'T' = 'A'
        update 'C' = 'G'
        update 'G' = 'C'

-- hamming distance between two strands
-- 'Error' between two DNA strands
hamm :: String -> String -> Int
hamm strand1 strand2 = foldr compareBases 0 zippedStrands
  where zippedStrands = zip strand1 strand2
        compareBases (base1, base2) errorCount | base1 == base2 = errorCount
                                               | base1 /= base2 = errorCount + 1

-- find substring indexes for repeats
-- Finding repeats in DNA strands
-- Strand -> Repeat -> List of indices of matches
-- Matching indices should be +1 in the final output (strings start at 1, not 0)
subs :: String -> String -> [Int]
subs strand repeat =  map (+1) $ filter repeatAtIndex possibleRepeatIndices 
  where repeatLength          = length repeat
        possibleRepeatIndices = findIndices (\char -> char == head repeat) strand
        repeatAtIndex index   = (take repeatLength $ drop index strand) == repeat

-- Convert an RNA (mRNA) string into a protein string
prot :: String -> String
prot =  concatMap rnaCodonTable . init . chunksOf 3

-- Computing GC content
--gc :: FullStringWithNewLinesButNo> -> (Name, DNAString)
-- [A, C, G, T]
gc :: String -> (String, String)
gc fullLine = (name, gcCount)
  where splitLine = splitOn "\n" fullLine
        name = head splitLine
        string = concat $ tail splitLine
        [a, c, g, t] = dna string
        gcCount = show $ 100 * ((fromIntegral (g + c)) / (fromIntegral (a + c + g + t)))

gcRun :: IO (String, String) 
gcRun = do
  contents <- readFile "rosalind_gc.txt"
  let fastaNameGCCounts = map gc $ tail $ splitOn ">" contents
  let maxFasta = maximumBy (\(_, x) (_, y) -> compare x y) fastaNameGCCounts
  putStrLn $ fst maxFasta
  putStrLn $ snd maxFasta
  return maxFasta 
        
-- Finding rabbit population after n months with k reproductive rate
-- 33 2 -> 2863311531
-- 5 3 -> 19
fib :: Int -> Int -> Int
fib 1 _ = 1
fib 2 _ = 1
fib months rate = (fib (months - 1) rate) + ((fib (months - 2) rate) * rate)

-- What is the probability of having a dominant allele
-- k => AA : m => Aa : n => aa
iprb :: Double -> Double -> Double -> Double
iprb k m n = sum [kk, km, kn, mk, mm, mn, nk, nm, nn]
  where t = k + m + n
        kk = (k / t) * ((k-1) / (t-1)) * domk
        km = (k / t) * (m     / (t-1)) * domk
        kn = (k / t) * (n     / (t-1)) * domk
        mk = (m / t) * (k     / (t-1)) * domk
        mm = (m / t) * ((m-1) / (t-1)) * dommm
        mn = (m / t) * (n     / (t-1)) * dommn
        nk = (n / t) * (k     / (t-1)) * domk
        nm = (n / t) * (m     / (t-1)) * dommn
        nn = (n / t) * ((n-1) / (t-1)) * domnn
        domk  = 1.0
        dommm = 0.75
        dommn = 0.5
        domnn = 0.0

rnaCodonTable :: String -> String
rnaCodonTable "UUU" = "F"    
rnaCodonTable "UUC" = "F"
rnaCodonTable "UUA" = "L"
rnaCodonTable "UUG" = "L"
rnaCodonTable "UCU" = "S"
rnaCodonTable "UCC" = "S"
rnaCodonTable "UCA" = "S"
rnaCodonTable "UCG" = "S"
rnaCodonTable "UAU" = "Y"
rnaCodonTable "UAC" = "Y"
rnaCodonTable "UAA" = "Stop"
rnaCodonTable "UAG" = "Stop"
rnaCodonTable "UGU" = "C"
rnaCodonTable "UGC" = "C"
rnaCodonTable "UGA" = "Stop"
rnaCodonTable "UGG" = "W"
rnaCodonTable "CUU" = "L"
rnaCodonTable "CUC" = "L"
rnaCodonTable "CUA" = "L"
rnaCodonTable "CUG" = "L"
rnaCodonTable "CCU" = "P"
rnaCodonTable "CCC" = "P"
rnaCodonTable "CCA" = "P"
rnaCodonTable "CCG" = "P"
rnaCodonTable "CAU" = "H"
rnaCodonTable "CAC" = "H"
rnaCodonTable "CAA" = "Q"
rnaCodonTable "CAG" = "Q"
rnaCodonTable "CGU" = "R"
rnaCodonTable "CGC" = "R"
rnaCodonTable "CGA" = "R"
rnaCodonTable "CGG" = "R"
rnaCodonTable "AUU" = "I"
rnaCodonTable "AUC" = "I"
rnaCodonTable "AUA" = "I"
rnaCodonTable "AUG" = "M"
rnaCodonTable "ACU" = "T"
rnaCodonTable "ACC" = "T"
rnaCodonTable "ACA" = "T"
rnaCodonTable "ACG" = "T"
rnaCodonTable "AAU" = "N"
rnaCodonTable "AAC" = "N"
rnaCodonTable "AAA" = "K"
rnaCodonTable "AAG" = "K"
rnaCodonTable "AGU" = "S"
rnaCodonTable "AGC" = "S"
rnaCodonTable "AGA" = "R"
rnaCodonTable "AGG" = "R"
rnaCodonTable "GUU" = "V"
rnaCodonTable "GUC" = "V"
rnaCodonTable "GUA" = "V"
rnaCodonTable "GUG" = "V"
rnaCodonTable "GCU" = "A"
rnaCodonTable "GCC" = "A"
rnaCodonTable "GCA" = "A"
rnaCodonTable "GCG" = "A"
rnaCodonTable "GAU" = "D"
rnaCodonTable "GAC" = "D"
rnaCodonTable "GAA" = "E"
rnaCodonTable "GAG" = "E"
rnaCodonTable "GGU" = "G"
rnaCodonTable "GGC" = "G"
rnaCodonTable "GGA" = "G"
rnaCodonTable "GGG" = "G"

