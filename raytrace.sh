#! /bin/bash
# create slice images for enriched b6 publication
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

usage ()
{
	echo "Usage: cat coords.txt | $(basename $0) -b BGFILE -l LABELS -s SURFACE -o OUTBASE  STAT1 THRESH1.."
	echo "  -b BGFILE           background "
	echo "  -l LABELS           labels  "
	echo "  -s SURFACE          surfaces, i.e. outline "
	echo "  -o OUTBASE          basename of output, coordinates get appended "
	echo "  STAT1 THRESH1       pairs of statistical overlay maps with associated (lower) threshold"
}
ycoord_to_png ()
{
	NPIXELS=600
	OUTXDIM=`python -c "print ${NPIXELS}.0/(5/3.0)"`
	OUTYDIM=`python -c "print ${NPIXELS}.0/2.0"`
	#echo output image dimensions: $OUTXDIM x $OUTYDIM

	ADD="-size $NPIXELS $NPIXELS -bg green -front -sup 1"   # sup > 1 gives green seems !!!
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
		ray_trace -output ${OBASE}_bg.rgb -under transparent -gray 100 1300 $BGFILE 2 1 ${OBASE}.obj $ADD

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
		convert -gravity Center -crop ${OUTXDIM}x${OUTYDIM}+0+0 +repage ${OBASE}_stat_on_rest.png ${OUTBASE}${SUFF}.png

	done
	echo
}


###############################################################################
# MAIN
###############################################################################

while getopts b:l:s:T:o:v opt
do
	case "$opt" in
		b)  BGFILE="$OPTARG";;
		l)  LABELS="$OPTARG";;
		s)  LABELSURFS="$OPTARG";;
		#i)  STATS="$OPTARG";;
		#t)  LOWERTHRESH="$OPTARG";;
		T)  UPPERTHRESH="$OPTARG";;
		o)  OUTBASE="$OPTARG";;
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


ARGS=($@)

TDIR=$TMPDIR/$$;
mkdir -p $TDIR
trap "{ rm -fr $TDIR; echo 'cleaned up temp. directory'; exit 255; }" SIGINT SIGTERM


while read YCOORD; do
	ycoord_to_png $YCOORD
done


rm -rf $TDIR
exit 0




# cd /projects/mice/jscholz/rot/stats_hr_new2+flip/dti_to_basket

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



# create cluster-extend thresholded signed images
# minccalc -expression 'A[0]>56?A[1]:0' clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005.mnc dti_meanCmeanFA_lmerp_mean-p-sign_to_basket.mnc clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005_gt56.mnc
# minccalc -expression 'A[0]>56?A[1]:0' clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005.mnc  dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket.mnc  clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005_gt56.mnc
# minccalc -expression 'A[0]>56?A[1]:0' clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1.mnc  final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign.mnc clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc
# minccalc -expression 'A[0]>56?A[1]:0' clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1.mnc  final_meanPsexRsideRid_lmerp_mean-p-corr-sign.mnc clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1_gt56.mnc 


#awk 'NR>1&&NR<5 {print $4}' clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005_gt56.mnc -t 0.995 -o slices/dti_group

#awk 'NR>1&&NR<44 {print $4}' clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005_gt56.mnc -t 0.995 -o slices/dti_behav_mean

#awk 'NR>1&&NR<11 {print $4}' clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -o slices/vol_beahv_mean clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1_gt56.mnc 0.9
######awk 'NR>1&&NR<11 {print $4}' clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1_gt56.mnc -t 0.9 -o slices/vol_beahv_mean

#awk 'NR>1&&NR<21 {print $4}' clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -o slices/vol_group clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc 0.9
######awk 'NR>1&&NR<21 {print $4}' clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc -t 0.9 -o slices/vol_group

# paste this into svg file:
# i=0; for f in `ls slices-src/dti_behav_mean_y*png | sort -r`; do x=$(($i%3*210-500)); y=$(($i/3*150)); echo "<image xlink:href=\"$f\" x=\"$x\" y=\"$y\" height=\"170.8\" width=\"223.3\" />"; let i+=1;  done`

# 3.40
# z = 3.14


# hippocampus hotspot of changes
# printf " -0.63\n-0.69\n" | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -o slices/hotspot clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc 0.9 clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005_gt56.mnc 0.995 clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005_gt56.mnc 0.995


