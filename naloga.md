Uporabimo oznake $n_{40}$ za število 40-letnikov v populaciji, $ifr_{40}$ naj
bo verjetnost, da 40-letnik umre ob okužbi z virusom in $p$ zaščita pred smrtjo,
ki jo zagotavlja cepivo. Potem je

- število umrlih necepljenih 40-letnikov: $`fu_{40} = 0.75~n_{40} \cdot ifr_{40}`$,
- število umrlih cepljenih 40-letnikov: $`fv_{40} = 0.25~n_{40} \cdot ifr_{40} \cdot (1 - p)`$ in
- število umrlih cepljenih 80-letnikov: $`fv_{80} = 0.5~n_{40} \cdot 120~ifr_{40} \cdot (1 - p)`$.

Če je med umrlimi cepljenimi in necepljenimi razmerje 1:1, velja:

```math
  \frac{fv_{40} + fv_{80}}{fu_{40}} = 4
  \frac{0.25 \cdot (1 - p) + 0.5 \cdot 120 \cdot (1 - p)}{0.75}
= (1 - p) \frac{241}{3} = 4
  p = 1 - \frac{12}{241} \approx 95\%
```
