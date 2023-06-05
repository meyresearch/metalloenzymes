i = createAtom   ZN  Zn2+  2
set i    element Zn
set i    position { 0 0 0 }
r = createResidue ZN
add r i
ZN = createUnit ZN
add ZN r
saveOff ZN ./ZN.lib
quit
