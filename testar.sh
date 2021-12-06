gcc -lm ep1.c

testes=( "1A" "1B" "1C" "1D" "1E" "1F" "1G"
  "2A" "2B" "2C" "2D" "2E" "2F" "2G" "2H" "2I" "2J" "2K"
  "3A" "3B" "3C" "3D" "3E" "3F" "3G" "3H"
  "4A" "4B" "4C" )

for i in "${testes[@]}" 
do
  echo 'Resultado do programa:'
  ./a.out < entradas/entrada"$i".txt
  echo ''
  echo 'Resultado esperado:'
  cat saidas/saida"$i".txt
  echo '---------------------------------'
done
