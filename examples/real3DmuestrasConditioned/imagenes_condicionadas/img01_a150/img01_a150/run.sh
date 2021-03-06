# 
# - Cargar pesos iniciales en weightout.dat (este es el unico archivo que se necesita pre-generar, se puede usar mags3d.c)
#
# - Se genera un variograma target, definido en genTargetVariogram, el cual crea el archivo targetvariogram.dat y targetimage.m
#   (internamente, crea una lattice del tamaño input, usando un kernel target hard-codeado, escribe la imagen convolucionada en disco
#    y aplica gam sobre esa imagen. Despues postprocesa el variograma resultado y lo deja en targetvariogram.dat)
# - Se lee weighout.dat y se crean los pesos iniciales, junto a su initialdistanceweights.dat (valores de los pesos ploteados desde el origen hacia afuera).
# - Se convoluciona una imagen random con los pesos iniciales y se calcula la funcion de costo.
# - Empieza iteracion. 
#


#rm initialimage*.dat currentimage*.dat

MAGS3D_HOME=/home/operedo/mags3d-from-git/mags3d

export OMP_NUM_THREADS=4
export OMP_SCHEDULE="static,32"

nx=40
ny=60
#nx=50
#ny=50
#nx=30
#ny=30
#nx=100
#ny=100
nz=12
#nx=40
#ny=60
#nz=12
xlo=0.0
ylo=0.0
zlo=0.0
h=10.0
#a=18.0
#a=90.0
a=150.0
#modefwd=1 
# 1: solo imprime weightout.dat
modeinv=1 
# 0: no usa genTargetVariogram, usa archivo ya generado a priori (var experimental de muestras.dat)
# 1: usa genTagetVariogram
usenscore=0

#iters=12000
#updatetemp=600
#tol=0.01
#initialtemp=0.01

iters=6000
report=100
updatetemp=300
tol=0.0001
initialtemp=0.0001

## synthetic data
#nx=128
#ny=128
#nz=1
#xlo=0.0
#ylo=0.0
#zlo=0.0
#h=1.0
#a=9.0
#modefwd=1 
## 1: solo imprime weightout.dat
#modeinv=1 
## 0: no usa genTargetVariogram, usa archivo ya generado a priori (var experimental de muestras.dat)
## 1: usa genTargetVariogram
#usenscore=0



### Opcional: generar targetdistanceweights.dat inicial, lanzando una simulacion fwd y escribiendo en disco los valores
##modefwd=2 
##echo "./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}"
##./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}
#

## Primero generar weightout.dat inicial, lanzando una simulacion fwd y escribiendo en disco los pesos iniciales
#modefwd=3 
##echo "./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}"
#$MAGS3D_HOME/bin/mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}
##

# Lanzar una simulacion fwd con datos condicionantes
modefwd=4 
#echo "./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}"

mu=1.054
echo "$MAGS3D_HOME/bin/mags3dConditioned.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd} ${mu}"
$MAGS3D_HOME/bin/mags3dConditioned.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd} ${mu} 2> salida.txt & 
#

## Segundo generar targetdistanceweights.dat inicial, lanzando una simulacion fwd y escribiendo en disco
#modefwd=2 
##echo "./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}"
#$MAGS3D_HOME/bin/mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}


#
## Luego se lanza la resolucion inversa
#echo "/usr/bin/time ./sa-invmags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modeinv} ${usenscore} 2> salida.err.log > salida.out.log" 
##/usr/bin/time ./sa-invmags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modeinv} ${usenscore} 2> salida.err.log > salida.out.log & 
#/usr/bin/time ./sa-invmags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modeinv} ${usenscore} ${iters} ${updatetemp} ${tol} ${initialtemp} 2> salida.err.log > salida.out.log & 
#/usr/bin/time $MAGS3D_HOME/bin/sa-cuda-invmags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modeinv} ${usenscore} gamv.par  2> salida.err.log > salida.out.log 

#/usr/bin/time $MAGS3D_HOME/bin/sa-cuda-invmags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modeinv} ${usenscore} ${iters} ${report} ${updatetemp} ${tol} ${initialtemp}  2> salida.err.log > salida.out.log 
#




#modefwd=3 
#echo "./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}"
#./mags3d.exe ${nx} ${ny} ${nz} ${xlo} ${ylo} ${zlo} ${h} ${a} ${modefwd}



