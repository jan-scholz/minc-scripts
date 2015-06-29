#! /bin/bash
# create slice images
# 2014-07-03,  jan.scholz@mouseimaging.ca
#
# TODO:
#  - make 'y' an argument


DIRECTION=y

LABELSCC=raytrace_colorcoding.txt 
POSCC=(hotred_new  cc_green cc_yellow cc_brown)
NEGCC=(hotblue_new cc_brown cc_brown)
LOWERTHRESH=0
UPPERTHRESH=1
NPIXELS=600

usage ()
{
	echo "Usage: cat coords.txt | $(basename $0) -b BGFILE -l LABELS -s SURFACE -o OUTBASE  STAT1 THRESH1.."
	echo "  -b BGFILE           background "
	echo "  -B BGTHRESH         threshold for background "
	echo "  -l LABELS           labels  "
	echo "  -s SURFACE          surfaces, i.e. outline "
	echo "  -o OUTBASE          basename of output, coordinates get appended "
	echo "  STAT1 THRESH1       pairs of statistical overlay maps with associated (lower) threshold"
}

coord_to_png ()
{
	#echo output image dimensions: $OUTXDIM x $OUTYDIM

    # removed -front because that doesn't work with direciton != y
	ADD="-size $NPIXELS $NPIXELS -bg green -sup 1"   # sup > 1 gives green seems !!!
	COORDS=`python -c "print ' '.join(['%.3f' % (${1}+i*0.056) for i in range(-1,2)])"`

	for COORD in $COORDS; do
		[ -z "$VERBOSE" ] || echo -n "$COORD "
		SUFF=`printf "_%s%+0.3f" $DIRECTION $COORD`
		OBASE=$TDIR/slice$SUFF
		make_slice $BGFILE ${OBASE}.obj $DIRECTION w $COORD
		set_object_colour ${OBASE}.obj ${OBASE}.obj green

		COMMAND="-output ${OBASE}.rgb -nolight"

		###############################################################################
		# SURFACES
		###############################################################################
		if [ ! -z "$LABELSURFS" ]; then
			plane_polygon_intersect $LABELSURFS $TDIR/lines.obj $DIRECTION $COORD
			set_object_colour $TDIR/lines.obj  $TDIR/lines.obj black
			ray_trace -output ${OBASE}_lines.rgb -nolight -line_width 0.03 $TDIR/lines.obj ${OBASE}.obj $ADD
		fi

		###############################################################################
		# BACKGROUND
		###############################################################################
		ray_trace -output ${OBASE}_bg.rgb -under transparent -gray $BGTHRESH $BGFILE 2 1 ${OBASE}.obj $ADD

		###############################################################################
		# LABELS
		###############################################################################
		if [ ! -z "$LABELS" ]; then
			ray_trace -output ${OBASE}_label.rgb -under transparent -green 0.05 1.8 $LABELS -1 1 ${OBASE}.obj $ADD
		fi

		###############################################################################
		# STATS
		###############################################################################
		for ((i=0;i<${#ARGS[@]}/3;i++)); do
			STATS=${ARGS[i*3]}
			LOWERTHRESH=${ARGS[i*3+1]}
			UPPERTHRESH=${ARGS[i*3+2]}
			
			ray_trace -output ${OBASE}_stat.rgb -under transparent -usercc ${POSCC[i]}  $LOWERTHRESH  $UPPERTHRESH $STATS 0 1 -under transparent -usercc ${NEGCC[i]} -$LOWERTHRESH -$UPPERTHRESH $STATS 0 1 ${OBASE}.obj $ADD

		done

		for f in ${OBASE}_*.rgb; do
			ALPHA=0%
			BC=0
			[ "${f##*_}" = "label.rgb" ] && ALPHA=-70%   # -80% almost fully transparent
			[ "${f##*_}" = "stat.rgb" ]  && ALPHA=-20%   # -80% almost fully transparent
			[ "${f##*_}" = "lines.rgb" ] && ALPHA=-30%
			[ "${f##*_}" = "bg.rgb" ]    && BC=0x0
			convert -rotate 180 -fill white -fuzz 10% +opaque '#0F0' -colorspace gray -threshold 90% -brightness-contrast $ALPHA $f ${f/.rgb/_mask.png}
			convert -rotate 180 -brightness-contrast $BC $f ${f/.rgb/.png}
		done

		convert -fuzz 10% -opaque '#0F0' -background black ${OBASE}_bg.png ${OBASE}_bg.png

		if [ ! -z "$LABELS" ]; then
			composite ${OBASE}_label.png ${OBASE}_bg.png ${OBASE}_label_mask.png ${OBASE}_bg.png
		fi
		
		if [ ! -z "$LABELSURFS" ]; then
			composite ${OBASE}_lines.png ${OBASE}_bg.png ${OBASE}_lines_mask.png ${OBASE}_bg.png
		fi
		
		composite ${OBASE}_stat.png ${OBASE}_bg.png ${OBASE}_stat_mask.png ${OBASE}_stat_on_rest.png


        ZDIM=`mincinfo -dimlength zspace $BGFILE`
        YDIM=`mincinfo -dimlength yspace $BGFILE`
        XDIM=`mincinfo -dimlength xspace $BGFILE`

        if [ $DIRECTION = x ]; then
            OUTXDIM=`python -c "print ${NPIXELS}.0*0.85"`
            OUTYDIM=`python -c "print ${NPIXELS}.0*${ZDIM}/${YDIM}*0.85"`
        fi
        if [ $DIRECTION = y ]; then
            OUTXDIM=`python -c "print ${NPIXELS}.0*0.85"`
            OUTYDIM=`python -c "print ${NPIXELS}.0*${ZDIM}/${XDIM}*0.85"`
        fi
        if [ $DIRECTION = z ]; then
            OUTXDIM=`python -c "print ${NPIXELS}.0*${XDIM}/${YDIM}*0.85"`
            OUTYDIM=`python -c "print ${NPIXELS}.0*0.85"`
        fi
        convert -gravity Center -crop ${OUTXDIM}x${OUTYDIM}+0+0 +repage ${OBASE}_stat_on_rest.png ${OUTBASE}${SUFF}.png
        
	done
	echo
}


###############################################################################
# MAIN
###############################################################################

while getopts b:B:l:s:T:o:d:v opt
do
	case "$opt" in
		b)  BGFILE="$OPTARG";;
		B)  BGTHRESH="$OPTARG";;
		l)  LABELS="$OPTARG";;
		s)  LABELSURFS="$OPTARG";;
		#i)  STATS="$OPTARG";;
		#t)  LOWERTHRESH="$OPTARG";;
		T)  UPPERTHRESH="$OPTARG";;
		o)  OUTBASE="$OPTARG";;
		d)  DIRECTION="$OPTARG";;
		v)  VERBOSE=1;;
		\?)  usage; exit 1;;
	esac
done
shift $(expr $OPTIND - 1)

[ -z "$TMPDIR" ] && { echo "ERROR: \$TMPDIR not set"; exit 1; }
[ "$( tty )" = 'not a tty' ] || { usage; exit 1; }

[ -f "$BGFILE" ]     || { echo "ERROR: could not find BACKGROUND_FILE: $BGFILE"; exit 1; }
[ -z "$LABELS" -o -f "$LABELS" ] || { echo "ERROR: could not find LABELS: $LABELS"; exit 1; }
[ -z "$LABELSURFS" -o -f "$LABELSURFS" ] || { echo "ERROR: could not find LABELSURFS: $LABELSURFS"; exit 1; }
#[ -f "$STATS" ]     || { echo "ERROR: could not find STATS: $STATS"; exit 1; }

for f in $LABELSCC ${POSCC[@]} ${NEGCC[@]}; do
	[ -f "$f" ] || { echo "ERROR: could not find file: $f"; exit 1; }
done


if [ -z "$BGTHRESH" ]; then
    P=(`mincstats -quiet -mean -stddev $BGFILE`)
    BGTHRESH=`python -c "print int(max(${P[0]} - 2*${P[1]},0)), int(max(${P[0]} + 2*${P[1]},1))"`
fi


ARGS=($@)

TDIR=$TMPDIR/$$;
mkdir -p $TDIR
trap "{ rm -fr $TDIR; echo 'cleaned up temp. directory'; exit 255; }" SIGINT SIGTERM


while read COORD; do
	coord_to_png $COORD
done

# echo $TDIR
# read

rm -rf $TDIR
exit 0


# Website with all the information:
# 
# http://www.bic.mni.mcgill.ca/~david/Ray_trace/ray_trace_tutorial.html

# Content of /tmp/colour_coding
# 
# 1.000 1.0000 0.0000 0.0000 
# 2 0 1 0
# 180  0 0 0
# 190 0 0 1
# 255 1 1 1
