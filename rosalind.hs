import Data.List

normalize :: [Int] -> String
normalize = intercalate " " . map show

-- [A, G, C, T]
-- Counts of nucleobases
dna :: String -> String
dna = normalize . foldr (\letter acc -> update letter acc) [0,0,0,0]
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
hamm s1 s2 = foldr (\(v1, v2) acc -> if v1 /= v2 then acc + 1 else acc) 0 (zip s1 s2)

-- find substring indexes for repeats
-- Finding repeats in DNA strands
-- Strand -> Repeat -> List of indices of matches
subs :: String -> String -> String
subs s r = normalize $ map (+1) $ filter repeatAtIndex possibleRepeatIndices 
  where rLength               = length r
        possibleRepeatIndices = findIndices (\v -> v == head r) s
        repeatAtIndex index   = (take rLength $ drop index s) == r 
