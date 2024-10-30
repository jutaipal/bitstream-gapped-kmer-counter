#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>

/* VERSION 0.1, 30 Oct 2024   part of AUTOSEED code modified to enable bitstream input (sequences of arbitrary length)   */
/* bitstream and multithreaded code written with assistance from Claude 3.5 Sonnet and ChatGPT4 and o1                   */
/* Debugged manually by J Taipale as bit shifting is not their forte                                                     */
/* Compile with clang -march=native -ffast-math -O3 -o bitstream_gkcount bitstream_gkcount.c                             */

/* GLOBAL VARIABLES */
uint64_t mask_ULL[42][42];   /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
uint64_t lowmask_ULL[42];    /* LOW MASK FOR EACH KMER */
uint64_t highmask_ULL[42];   /* HIGH MASK FOR EACH KMER */
short int Nlength = 32;              /* LENGTH OF MASK */

#define NUM_FILES 2
#define MAX_GAP_LENGTH 10
#define KMER_VARIATION 1
#define BUFFER_SIZE 1024
#define MAX_KMER_LEN 32  // Adjust if needed for your use case
#define MAX_SEQ_LEN 128  // Adjust if needed for your use case



/* SUBROUTINE THAT DETERMINES IF GAP IS AT EITHER OF THE CENTER POSITIONS, IF count_also is set to != 1 returns true */
short int Centergap (short int count_also_spaced_kmers, short int kmer_length, short int gap_position)
{
if (count_also_spaced_kmers != 1) return (1);
if (gap_position == kmer_length / 2) return (1);
if (kmer_length % 2 == 1 && gap_position == kmer_length / 2 - 1) return (1);
else return (0);
}

/* SUBROUTINE THAT DETERMINES IF A GAPPED KMER IS A LOCAL MAXIMUM WITHIN HUDDINGE DISTANCE OF 1 */
/* SEE NITTA ET AL. eLIFE 2015 Methods and Supplementary Figure 1 for algorithm description     */
/* https://doi.org/10.7554/eLife.04837.004                                                      */
/* Bases encoded as bits, A = 00, C = 01, G = 10, T = 11                                        */

short int Localmax(long int *****results, short int file_number, short int current_kmer_length, short int shortest_kmer, short int longest_kmer_counted, short int current_gap_position, short int current_gap_length, long int current_kmer, double kmer_length_difference_cutoff)
{
    short int count_also_spaced_kmers = 1;
    short int too_long_kmer = longest_kmer_counted + 1;

    long int kmer1_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
    long int kmer2_incidence;
    long int compared_kmer = current_kmer;
    signed short int position;
    short int counter;
    short int first_half;
    long int lowbit = 1;  // "01" used to toggle low bit of base encoding (between A-C or G-T), shifted to correct base position in loops
    long int highbit = 2; // "10" used to toggle high bit of base encoding (between A-G or C-T), shifted to correct base position in loops
    short int shift;
    short int true_gap_position = current_gap_position;
    short int true_gap_length = current_gap_length;
    short int start;
    short int end;
    short int left = 0;
    short int right = 1;
    signed short int position2;
    
    /* Substitution; HAMMING OF 1, RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT */
    for(position=0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2)
    {
        compared_kmer = lowbit ^ current_kmer;  // SHIFTS between A-C or G-T of original kmer at position
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);  // IF SUBSTITUTED KMER HAS HIGHER COUNTS, RETURNS 0 TO INDICATE THAT THE QUERY KMER IS NOT LOCAL MAXIMUM
        compared_kmer = highbit ^ current_kmer; // SHIFTS between A-G or C-T of original kmer at position
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer; // SHIFTS between A-C or G-T of previously shifted kmer to cover last of possible 3 substitutions (total effect is A-T, C-G)
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    /* Shift; FULL SHIFT BY ONE */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    /* ONLY LOOK AT FULL SHIFT FOR UNGAPPED KMERS, OR FOR GAPPED KMERS IF GAP CAN SHIFT ALSO (ALL GAP POSITIONS HAVE BEEN COUNTED) */
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position < current_kmer_length) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1))
    {
    /* COMPARED FULL SHIFT RIGHT (same shift in gap if any), RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT  */
    compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit >>= 2; highbit >>= 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position > 0) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))
    {
    /* COMPARED FULL SHIFT LEFT (same shift in gap if any) */
    compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit = 1; highbit = 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position += shift;
    lowbit = 1; highbit = 2;
    if (count_also_spaced_kmers != 0)
    {
    if(current_gap_position == 0) current_gap_position = current_kmer_length / 2;
    
    /* Longer Gap; COMPARE TO KMER WITH LONGER GAP */
    current_gap_length++;
    if (current_gap_length < Nlength - current_kmer_length && current_gap_length < MAX_GAP_LENGTH && true_gap_length != 0)
    {
    compared_kmer = (current_kmer & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    
    /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
    for(first_half = 0; first_half < 2; first_half++)
    {
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    
    /* SWITCHES TO FIRST HALF */
    compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
    }
    }
        
    /* Shorter Gap; COMPARE TO KMER WITH SHORTER GAP */
    current_gap_length--;
    current_gap_length--;

    lowbit = 1; highbit = 2;
    lowbit <<= ((current_kmer_length - true_gap_position - 1)*2); highbit <<= ((current_kmer_length - true_gap_position - 1)*2);
    if (current_gap_length == 0) current_gap_position = 0;
    if (current_gap_length >= 0)
    {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer >> 2) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
        /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
        for(first_half = 0; first_half < 2; first_half++)
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            
            /* SWITCHES TO FIRST HALF */
            compared_kmer = lowmask_ULL[current_kmer_length-1] & ((current_kmer << 2) & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
            lowbit <<= 2; highbit <<= 2;
        }
    }
        
        /* COMPARES DIFFERENT GAP POSITIONS (SHIFTED BY ONE) */
        current_gap_length = true_gap_length;
        if ((count_also_spaced_kmers == 2 || (count_also_spaced_kmers == 1 && current_kmer_length % 2 == 1)) && true_gap_length > 0)
        {
        current_gap_position = true_gap_position + 1;
        lowbit = 1; highbit = 2;
        lowbit <<= ((current_kmer_length - current_gap_position)*2); highbit <<= ((current_kmer_length - current_gap_position)*2);
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));

        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH RIGHT BY ONE */
        for(first_half = 0; first_half < 2; first_half++)
        {
        if (current_gap_position < current_kmer_length && current_gap_position > 0 && (count_also_spaced_kmers == 2 ||
        (count_also_spaced_kmers == 1 && ((current_gap_position == current_kmer_length / 2) || current_gap_position == current_kmer_length / 2 + 1))))
        {
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        }
            
        /* SWITCHES LEFT BY ONE */
        current_gap_position--;
        current_gap_position--;
        compared_kmer = ((current_kmer) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
        lowbit <<= 2; highbit <<= 2;
        }
        }
        
        start = 1;
        end = current_kmer_length;
        
        /* COMPARES UNGAPPED KMER TO ALL SINGLE GAPS IN ALL POSITIONS */
        if (count_also_spaced_kmers != 0 && true_gap_length == 0)
        {
            if(count_also_spaced_kmers == 1)
            {
                start = current_kmer_length / 2;
                end = start + 1 + (current_kmer_length % 2);
            }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
            for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
            {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
           
        /* END */
        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
        }
        }
        }
        
    }
/* Shorter; COMPARES KMER WITH ONE SHORTER */
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
current_kmer_length--;
end = current_kmer_length;
if (current_kmer_length >= shortest_kmer)
{
/* IF NO GAP, INSERTS GAP AT ALL ALLOWED POSITIONS */
if(current_gap_length == 0)
{
    if (count_also_spaced_kmers != 0)
    {
        if(count_also_spaced_kmers == 1)
        {
            start = current_kmer_length / 2;
            end = start + 1 + (current_kmer_length % 2);
        }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                if(kmer2_incidence  * kmer_length_difference_cutoff > kmer1_incidence) return (0);
        }
    }
}

current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
if (count_also_spaced_kmers != 1) {left = 1; right = 1;}
if (current_gap_position == current_kmer_length / 2 && current_kmer_length % 2 == 0 && count_also_spaced_kmers == 1) {left = 1; right = 0;}
if (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1) left = 1;
    
/* LEFT PART */
if (current_gap_position < current_kmer_length)
{
    if (left == 1 || true_gap_length == 0)
    {
compared_kmer = (current_kmer >> 2);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    }
/* RIGHT PART */
if (current_gap_position != 1 && right == 1)
{
if (current_gap_position > 0) current_gap_position--;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_position = true_gap_position;
/* Shorter with Longer Gap; LONGER GAP */
if(count_also_spaced_kmers != 0 && true_gap_position != 0 && current_gap_length < MAX_GAP_LENGTH)
{
current_gap_length++;
if (current_gap_position < current_kmer_length && left == 1)
{
/* RIGHT BASE GAPPED */
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    /* printf("DEBUG Localmax: k=%d, gap_pos=%d, gap_len=%d, indices: high=%d, low=%d\n",
           current_kmer_length, current_gap_position, current_gap_length,
           current_kmer_length - current_gap_position,
           current_kmer_length-current_gap_position-1); fflush(stdout); */
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
/* LEFT BASE GAPPED */
if(current_gap_position > 1 && right == 1)
{
current_gap_position--;
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_length--;
current_gap_position++;
}
    /* COMPARES HANGING SINGLE BASE TO UNGAPPED KMER */
    current_gap_position = true_gap_position;
    current_gap_length = true_gap_length;
    if(current_gap_position == 1)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    else if(current_kmer_length-current_gap_position == 0)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer >> 2);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
/* Longer; COMPARES KMER WITH ONE LONGER */
current_kmer_length++;
current_kmer_length++;
if (current_kmer_length < too_long_kmer)
{
compared_kmer = (current_kmer << 2);
/* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
{
if(count_also_spaced_kmers != 1 || true_gap_length == 0 || (first_half == 1 || (count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)))
{
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
   // printf("compared kmer count: %li\n", kmer2_incidence);
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = lowbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = highbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = lowbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
}
/* END */
compared_kmer = current_kmer;
lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
if (true_gap_length == 0) continue;
current_gap_position++;
if (current_gap_position > current_kmer_length || (count_also_spaced_kmers == 1 && current_gap_position != (current_kmer_length / 2 + (current_kmer_length % 2)))) break;
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    

    
/* Longer with Shorter Gap; COMPARES TO ONE SHORTER GAP LENGTH */
if(count_also_spaced_kmers != 0 && true_gap_length >= 1)
{
current_gap_length--;
lowbit = 1; highbit = 2;
lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
/* NO GAP LEFT */
if (current_gap_length == 0) current_gap_position = 0;
    
    /* LOOP TO ANALYZE ADDED BASE ON EITHER SIDE OF GAP, STARTS WITH RIGHT */
    for(first_half = 0; first_half < 2; first_half++)
    {
        if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position)))
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
        }
        
        /* SWITCHES ADDED BASE TO LEFT */
        current_gap_position++;
        if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
        compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    }
current_gap_length++;
}
}
current_gap_position = true_gap_position;
return (1);
}

/* SUBROUTINE THAT GENERATES 128 bit BITMASKS FOR NUCLEOTIDES, KMER STRINGS AND DELETIONS */
void GenerateMask ()
{
    short int counter;
    short int start_position;
    short int position;
    short int current_kmer_length;
    /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
    for (mask_ULL[1][0] = 3ULL, counter = 0; counter < Nlength; counter++) mask_ULL[1][counter+1] = mask_ULL[1][counter] << 2;
    /* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING */
    for(current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++)
    {
        for (start_position = 0; start_position < Nlength-current_kmer_length; start_position++)
        {
            for (mask_ULL[current_kmer_length][start_position]=mask_ULL[1][start_position], position = start_position+1; position < current_kmer_length+start_position; position++) {mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position];}
        }
    }
    /* GENERATES HIGH AND LOW mask_ULLS FOR DELETIONS */
    for (lowmask_ULL[0] = mask_ULL[1][0], position = 1; position < Nlength-2; position++) lowmask_ULL[position] = lowmask_ULL[position-1]+mask_ULL[1][position];
    for (highmask_ULL[Nlength-2] = mask_ULL[1][Nlength-2], position = Nlength-3; position > 0; position--) highmask_ULL[position] = highmask_ULL[position+1]+mask_ULL[1][position];
}

// LineBuffer structure to include 5D results
typedef struct {
    uint64_t***** results;  // 5D array [file][kmer_len][gap_pos][gap_len][kmer]
    uint64_t* bitstream;
    int base_kmer_length;
    int max_seq_length;
    int file_index;
    char finished;
    pthread_mutex_t mutex;
    pthread_cond_t not_full;
    pthread_cond_t not_empty;
    char** sequences;
    size_t* seq_lengths;
    int write_pos;
    int read_pos;
    int count;
    int shortest_kmer;     // Add these
    int longest_kmer;      // two fields
    FILE* file;
} LineBuffer;

// Function declarations for externally defined functions
// Function to encode a DNA character to 2-bit representation.
uint64_t encode_base(char base) {
    switch (base) {
        case 'A': return 0x0ULL; // 00
        case 'C': return 0x1ULL; // 01
        case 'G': return 0x2ULL; // 10
        case 'T': return 0x3ULL; // 11
        default: {/*printf("\nError in DNA sequence");*/return UINT64_MAX;} // Return max if invalid character
    }
}

// Function to convert a kmer index back to a DNA string.
void decode_kmer(uint64_t kmer, size_t kmer_length, size_t gap_position, size_t gap_length, char *out) {
    int Ncount = gap_length;
    int kmer_position = kmer_length;
    int i = kmer_position + gap_length - 1;

    for (; i >= 0; i--) {
        if (kmer_position == gap_position && Ncount > 0) {
            out[i] = 'n';
            Ncount--;
        }
        else {
            switch (kmer & 0x3) {
                case 0x0: out[i] = 'A'; break;
                case 0x1: out[i] = 'C'; break;
                case 0x2: out[i] = 'G'; break;
                case 0x3: out[i] = 'T'; break;
            }
            kmer >>= 2;
            kmer_position--;
        }
    }
    out[kmer_length+gap_length] = '\0'; // Null-terminate the string.
}

void print_bitstream_as_dna(const uint64_t *bitstream, size_t block_index, size_t dna_length) {
    printf("Block %zu: ", block_index);
    for (int i = 31; i >= 0; i--) {
        uint64_t mask = 0x3ULL << (2 * i);
        char base = ((32 - i) < dna_length) ?
                    "ACGT"[(bitstream[block_index] & mask) >> (2 * i)] :
                    'a';
        printf("%c", base);
    }
    printf("\n");
}


short int encode_dna_sequence(const char *dna_sequence, uint64_t *bitstream, size_t *bitstream_length) {
    size_t dna_length = strlen(dna_sequence);
    *bitstream_length = (2 * dna_length + 63) / 64;
    memset(bitstream, 0, (*bitstream_length) * sizeof(uint64_t));
    
    // printf("DEBUG: Encoding sequence of length %zu\n", dna_length);
    
    // Initial encoding
    size_t block;
    for (size_t i = 0; i < dna_length; i++) {
        uint64_t encoded_base = encode_base(dna_sequence[i]);
        if (encoded_base == UINT64_MAX) return (0); // DEFECTIVE SEQUENCE
        block = (*bitstream_length) - 1 - i * 2 / 64;
        size_t offset = 62 - (i * 2 % 64);
        bitstream[block] |= (encoded_base << offset);
    }
    bitstream[(*bitstream_length)] = 0;  // add 0 block after so that shifting last block will not bring in garbage seq

    // printf("\nAfter initial encoding:\n");
    // for (size_t i = 0; i < *bitstream_length; i++) {
        // print_bitstream_as_dna(bitstream, i, dna_length);
    // }

    // Handle short sequences
    if (dna_length <= 32) {
        bitstream[0] = bitstream[0] >> (64 - 2 * dna_length);
        // printf("\nAfter short sequence adjustment:\n");
        // print_bitstream_as_dna(bitstream, 0, dna_length);
    }
    
    // Handle longer sequences
    if (dna_length > 32 && dna_length % 32 != 0) {
        size_t left_shift = (dna_length * 2) % 64;
        size_t right_shift = 64 - left_shift;
        
        /*printf("\nBefore shifting blocks:\n");
        for (size_t i = 0; i < *bitstream_length; i++) {
            print_bitstream_as_dna(bitstream, i, dna_length);
            }*/
        
        uint64_t low_bits;
        //Shift blocks maintaining proper bit transitions
        for (size_t block = 0; block < (*bitstream_length); block++) {
            // Get high bits from next block that need to be brought into this block
            uint64_t high_bits = bitstream[block + 1] << left_shift;
            low_bits = (bitstream[block] >> right_shift);
            // Shift current block right and OR in the high bits from next block
            bitstream[block] = low_bits | high_bits;
        }
        
        /*printf("\nAfter shifting blocks:\n");
        for (size_t i = 0; i < *bitstream_length; i++) {
            print_bitstream_as_dna(bitstream, i, dna_length);
        }*/
        
    }
    return(1);
}

void count_kmers(uint64_t *results, uint64_t *bitstream, size_t sequence_length, size_t kmer_length, size_t gap_position, size_t gap_length) {
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    size_t total_kmers = sequence_length - kmer_length + 1 - gap_length;
    size_t counted_kmers = 0;
    size_t incoming_bits;
    size_t last_source_block_moved = 0;
    uint64_t kmer;
    size_t current_block = 0;
    // char* kmer_str = (char*)malloc(kmer_length + MAX_GAP_LENGTH + 1);
    // kmer_str[kmer_length] = '\0';
    
    // printf("DEBUG: Starting kmer counting, expecting %zu kmers\n", total_kmers);
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;

    uint64_t *current_bitstream = bitstream;

    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
        // printf("DEBUG: Loaded block %zu to main: %016lx\n", current_block++, main_reg);
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
        // printf("DEBUG: Loaded block %zu to source: %016lx\n", current_block++, source_reg);
    }

    size_t total_blocks = (sequence_length + 31) / 32;  // How many blocks we should have
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64; // counts bits left in source and finally main reg

//asm("kmer_count_loop_start:");
    while (counted_kmers < total_kmers) {
        // printf("\nDEBUG: --- Iteration %zu ---\n", counted_kmers);
        // printf("main_reg: %016lx, source_reg: %016lx\n", main_reg, source_reg);
        // printf("remaining_bits_in_main: %zu\n", remaining_bits_in_main);

        // Extract and count kmer
        if (gap_length == 0) {
            kmer = main_reg & ((1ULL << (2 * kmer_length)) - 1);
        } else {
            uint64_t mask = (1ULL << (2 * (kmer_length - gap_position))) - 1;
            kmer = ((main_reg & mask) ^
                   ((main_reg >> (2 * gap_length)) & ~mask)) &
                   ((1ULL << (2 * kmer_length)) - 1);
        }

        // decode_kmer(kmer, kmer_length, gap_position, gap_length, kmer_str);
        // printf("DEBUG: Kmer %zu = %s from block %zu\n", counted_kmers, kmer_str, current_block - 1);
                             
        results[kmer]++;
        counted_kmers++;

        // Shift bits
        uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
        // printf("DEBUG: new_bits: %016lx\n", new_bits);

        #ifdef __x86_64__
        asm volatile (
            "shrdq $2, %1, %0"
            : "+r" (main_reg)
            : "r" (source_reg)
            : "cc"
        );
        #else
        main_reg = (main_reg >> shift_amount) | new_bits;
        #endif

        source_reg >>= shift_amount;
        remaining_bits_in_main -= 2;

        // printf("DEBUG: After shift - main_reg: %016lx, source_reg: %016lx, blocks remaining:%li\n", main_reg, source_reg, total_number_of_blocks_remaining);

        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 && total_number_of_blocks_left_to_move_to_source == 0) {
            remaining_bits_in_main = (sequence_length % 32) * 2;
            last_source_block_moved = 1;
            // printf("DEBUG: last block in source. Remaining blocks, bits %li, %li\n", total_number_of_blocks_left_to_move_to_source, remaining_bits_in_main);
        }
        
        // Critical part: loading next block
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 && total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
                current_block++;
                // printf("DEBUG: new block moved to source. Remaining blocks %li, %li\n", total_number_of_blocks_left_to_move_to_source, remaining_bits_in_main);
            }

    }
//asm("kmer_count_loop_end:");
    
    // free(kmer_str);
    // printf("DEBUG: Counted %zu kmers\n", counted_kmers);
}

uint64_t***** allocate_5d_results(int base_kmer_len) {
    uint64_t***** results = malloc(NUM_FILES * sizeof(uint64_t****));
    
    for (int f = 0; f < NUM_FILES; f++) {
        // We need arrays only for k-1, k, and k+1
        results[f] = malloc((base_kmer_len + 2) * sizeof(uint64_t***));
        
        for (int k = base_kmer_len - 1; k <= base_kmer_len + 1; k++) {
            if (k <= 0) continue;
            
            // Calculate exact size needed for this kmer length
            uint64_t total_kmers = 1ULL << (2 * k);
            //printf("DEBUG: Allocating for k=%d, total_kmers=%lu\n", k, total_kmers);
            
            // Allocate only position 0 and center positions
            int center1, center2;
            if (k % 2 == 0) {
                center1 = k / 2;
                center2 = -1;  // No second center for even length
            } else {
                center1 = (k - 1) / 2;
                center2 = center1 + 1;
            }
            
            // Allocate array for gap positions, only need indices 0 and center position(s)
            results[f][k] = malloc((k + 1) * sizeof(uint64_t**));
            for (int g = 0; g <= k; g++) {
                results[f][k][g] = NULL;  // Initialize all to NULL
            }
            
            // Allocate position 0 for ungapped kmers
            results[f][k][0] = malloc(sizeof(uint64_t*));  // Only need index 0
            results[f][k][0][0] = calloc(total_kmers, sizeof(uint64_t));
            if (!results[f][k][0][0]) {
                fprintf(stderr, "Failed to allocate memory for ungapped kmers at k=%d (size=%llu)\n",
                        k, total_kmers);
                exit(1);
            }
            
            // Allocate first center position
            results[f][k][center1] = malloc((MAX_GAP_LENGTH + 1) * sizeof(uint64_t*));
            for (int l = 1; l <= MAX_GAP_LENGTH; l++) {  // Start from 1 for gaps
                results[f][k][center1][l] = calloc(total_kmers, sizeof(uint64_t));
                if (!results[f][k][center1][l]) {
                    fprintf(stderr, "Failed to allocate memory for gapped kmers at k=%d, center=%d, gap=%d (size=%llu)\n",
                            k, center1, l, total_kmers);
                    exit(1);
                }
            }
            
            // Allocate second center position if needed (odd length)
            if (center2 != -1) {
                results[f][k][center2] = malloc((MAX_GAP_LENGTH + 1) * sizeof(uint64_t*));
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {  // Start from 1 for gaps
                    results[f][k][center2][l] = calloc(total_kmers, sizeof(uint64_t));
                    if (!results[f][k][center2][l]) {
                        fprintf(stderr, "Failed to allocate memory for gapped kmers at k=%d, center=%d, gap=%d (size=%llu)\n",
                                k, center2, l, total_kmers);
                        exit(1);
                    }
                }
            }
        }
    }
    return results;
}

// Matching free function that only frees allocated memory
void free_5d_results(uint64_t***** results, int base_kmer_len) {
    if (!results) return;
    
    for (int f = 0; f < NUM_FILES; f++) {
        if (!results[f]) continue;
        
        for (int k = base_kmer_len - 1; k <= base_kmer_len + 1; k++) {
            if (k <= 0 || !results[f][k]) continue;
            
            // Free ungapped kmer storage
            if (results[f][k][0]) {
                free(results[f][k][0][0]);
                free(results[f][k][0]);
            }
            
            // Calculate center positions
            int center1, center2;
            if (k % 2 == 0) {
                center1 = k / 2;
                center2 = -1;
            } else {
                center1 = (k - 1) / 2;
                center2 = center1 + 1;
            }
            
            // Free first center position
            if (results[f][k][center1]) {
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                    free(results[f][k][center1][l]);
                }
                free(results[f][k][center1]);
            }
            
            // Free second center position if it exists
            if (center2 != -1 && results[f][k][center2]) {
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                    free(results[f][k][center2][l]);
                }
                free(results[f][k][center2]);
            }
            
            free(results[f][k]);
        }
        free(results[f]);
    }
    free(results);
}

void init_buffer(LineBuffer* buffer, int base_kmer_length, int file_index, const char* filename) {
    buffer->base_kmer_length = base_kmer_length;
    buffer->file_index = file_index;
    buffer->finished = 0;
    buffer->write_pos = 0;
    buffer->read_pos = 0;
    buffer->count = 0;
    buffer->shortest_kmer = buffer->base_kmer_length - 1;
    buffer->longest_kmer = buffer->base_kmer_length + 1;
    
    pthread_mutex_init(&buffer->mutex, NULL);
    pthread_cond_init(&buffer->not_full, NULL);
    pthread_cond_init(&buffer->not_empty, NULL);
    
    buffer->sequences = malloc(BUFFER_SIZE * sizeof(char*));
    buffer->seq_lengths = malloc(BUFFER_SIZE * sizeof(size_t));
    for (int i = 0; i < BUFFER_SIZE; i++) {
        buffer->sequences[i] = malloc(1024 * sizeof(char));  // Adjust size as needed
    }
    
    buffer->file = fopen(filename, "r");
    if (!buffer->file) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
}

void destroy_buffer(LineBuffer* buffer) {
    for (int i = 0; i < BUFFER_SIZE; i++) {
        free(buffer->sequences[i]);
    }
    free(buffer->sequences);
    free(buffer->seq_lengths);
    
    pthread_mutex_destroy(&buffer->mutex);
    pthread_cond_destroy(&buffer->not_full);
    pthread_cond_destroy(&buffer->not_empty);
    
    fclose(buffer->file);
}

/* Modified producer to handle empty lines
void* producer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    char line[1024];
    
    //printf("Producer: Starting\n");
    while (fgets(line, sizeof(line), buffer->file)) {
        // Skip empty lines or lines with just whitespace
        size_t len = strlen(line);
        if (len == 0) continue;
        
        // Trim trailing whitespace and newlines
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ')) {
            line[--len] = '\0';
        }
        
        // Skip if line is empty after trimming
        if (len == 0) continue;
        
        //printf("Producer: Got line of length %zu\n", len);
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == BUFFER_SIZE) {
            //printf("Producer: Buffer full (count=%d), waiting\n", buffer->count);
            pthread_cond_wait(&buffer->not_full, &buffer->mutex);
        }
        
        strcpy(buffer->sequences[buffer->write_pos], line);
        buffer->seq_lengths[buffer->write_pos] = len;
        buffer->write_pos = (buffer->write_pos + 1) % BUFFER_SIZE;
        buffer->count++;
        
        //printf("Producer: Added line at pos %d, count now %d\n",
        //       (buffer->write_pos - 1 + BUFFER_SIZE) % BUFFER_SIZE, buffer->count);
        
        pthread_cond_signal(&buffer->not_empty);
        pthread_mutex_unlock(&buffer->mutex);
    }
    
    //printf("Producer: Reached end of file\n");
    pthread_mutex_lock(&buffer->mutex);
    buffer->finished = 1;
    pthread_cond_broadcast(&buffer->not_empty);
    pthread_mutex_unlock(&buffer->mutex);
    //printf("Producer: Finished\n");
    
    return NULL;
} */

void* producer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    char line[1024];
    time_t last_print = time(NULL);
    
    while (fgets(line, sizeof(line), buffer->file)) {
        // Skip empty lines or lines with just whitespace
        size_t len = strlen(line);
        if (len == 0) continue;
        
        // Trim trailing whitespace and newlines
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ')) {
            line[--len] = '\0';
        }
        
        // Skip if line is empty after trimming
        if (len == 0) continue;
        
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == BUFFER_SIZE) {
            time_t current = time(NULL);
            if (current - last_print >= 10) {
                // printf("Producer: Still waiting for buffer space after %ld seconds. Buffer count=%d\n",
                //       current - last_print, buffer->count);
                // fflush(stdout);
                last_print = current;
            }
            pthread_cond_wait(&buffer->not_full, &buffer->mutex);
        }
        
        strcpy(buffer->sequences[buffer->write_pos], line);
        buffer->seq_lengths[buffer->write_pos] = len;
        buffer->write_pos = (buffer->write_pos + 1) % BUFFER_SIZE;
        buffer->count++;
        
        pthread_cond_signal(&buffer->not_empty);
        pthread_mutex_unlock(&buffer->mutex);
    }
    
    pthread_mutex_lock(&buffer->mutex);
    buffer->finished = 1;
    pthread_cond_broadcast(&buffer->not_empty);
    pthread_mutex_unlock(&buffer->mutex);
    
    return NULL;
}

void* consumer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    size_t line_count = 1;
    
    // printf("Consumer: Starting\n");
    while (1) {
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == 0) {
            if (buffer->finished) {
                pthread_mutex_unlock(&buffer->mutex);
                return NULL;
            }
            pthread_cond_wait(&buffer->not_empty, &buffer->mutex);
        }
        
        char* sequence = buffer->sequences[buffer->read_pos];
        size_t seq_length = buffer->seq_lengths[buffer->read_pos];
        
        buffer->read_pos = (buffer->read_pos + 1) % BUFFER_SIZE;
        buffer->count--;
        
        pthread_cond_signal(&buffer->not_full);

        // Keep mutex locked during sequence processing
        size_t bitstream_length;
        memset(buffer->bitstream, 0, 64 * sizeof(uint64_t));
        
        if (encode_dna_sequence(sequence, buffer->bitstream, &bitstream_length) == 0)
        {
            printf("\nError in DNA sequence on line %li, sequence rejected", line_count);
            pthread_mutex_unlock(&buffer->mutex);
            line_count++;
            continue;
        }
        
        pthread_mutex_unlock(&buffer->mutex);
        line_count++;

        // In the kmer counting part:
        // Process ungapped kmers within valid range
        for (int k = buffer->shortest_kmer; k <= buffer->longest_kmer; k++) {
            if (k <= 0) continue;
            if (buffer->results[buffer->file_index][k][0] &&
                buffer->results[buffer->file_index][k][0][0]) {
                count_kmers(buffer->results[buffer->file_index][k][0][0],
                           buffer->bitstream, seq_length, k, 0, 0);
            }
        }

        // Process kmers within valid gap position and gap length range
        for (int k = buffer->shortest_kmer; k <= buffer->longest_kmer; k++) {
            if (k <= 0) continue;
            
            // printf("DEBUG: Processing k=%d\n", k);
            
            // Calculate gap positions based on kmer length
            int number_of_gap_positions;
            int gap_positions[2];
            
            if (k % 2 == 0) {
                number_of_gap_positions = 1;
                gap_positions[0] = k/2;
                // printf("DEBUG: Even length kmer, gap position at %d\n", gap_positions[0]);
            } else {
                number_of_gap_positions = 2;
                gap_positions[0] = k/2;
                gap_positions[1] = k/2 + 1;
                // printf("DEBUG: Odd length kmer, gap positions at %d and %d\n",
                //        gap_positions[0], gap_positions[1]);
            }
            
            // For each gap position
            for (int i = 0; i < number_of_gap_positions; i++) {
                int gap_pos = gap_positions[i];
                // printf("DEBUG: Processing gap position %d\n", gap_pos);
                
                // For each gap length (1-10)
                for (int gap_len = 1; gap_len <= 10; gap_len++) {
                    // Skip if kmer + gap would be longer than sequence
                    if (k + gap_len > seq_length) {
                        // printf("DEBUG: Skipping k=%d gap_len=%d (total %d > seq_length %zu)\n",
                        //        k, gap_len, k + gap_len, seq_length);
                        continue;
                    }

                    // printf("DEBUG: Processing gap length %d\n", gap_len);
                    
                    if (buffer->results[buffer->file_index][k][gap_pos] &&
                        buffer->results[buffer->file_index][k][gap_pos][gap_len]) {
                        
                        // uint64_t array_size = 1ULL << (2 * k);
                        // printf("DEBUG: Result array size for k=%d: %" PRIu64 " entries (%" PRIu64 " bytes)\n",
                        //        k, array_size, array_size * sizeof(uint64_t));
                        
                        count_kmers(buffer->results[buffer->file_index][k][gap_pos][gap_len],
                                  buffer->bitstream, seq_length, k, gap_pos, gap_len);
                        
                        // printf("DEBUG: Finished counting for k=%d gap_pos=%d gap_len=%d\n",
                        //        k, gap_pos, gap_len);
                    }
                }
            }
        }
        
        // printf("Consumer: Finished processing sequence\n\n");
    }
}

// And modify print function to take file parameter:
void print_kmers_single_config(uint64_t***** results, int file, int kmer_len,
                             int shortest_kmer, int longest_kmer,
                             int gap_pos, int gap_len, double threshold,
                             int only_local_max, double length_diff_cutoff) {
    char* kmer_str = malloc(kmer_len + MAX_GAP_LENGTH + 1);
    if (!kmer_str) return;
    
    uint64_t total_kmers = 1ULL << (2 * kmer_len);
    
    // printf("\nResults for k=%d gap_pos=%d gap_len=%d:\n",
    //       kmer_len, gap_pos, gap_len);
    
    for (uint64_t kmer = 0; kmer < total_kmers; kmer++) {
        uint64_t signal_count = results[1][kmer_len][gap_pos][gap_len][kmer];
        uint64_t background_count = results[0][kmer_len][gap_pos][gap_len][kmer];

        if (signal_count >= threshold) {
            short int is_localmax = Localmax((long int*****)results,
                                           1, kmer_len,
                                           shortest_kmer, longest_kmer,
                                           gap_pos, gap_len, kmer,
                                           length_diff_cutoff);
            decode_kmer(kmer, kmer_len, gap_pos, gap_len, kmer_str);
            if (!only_local_max || is_localmax == 1) printf("%s\t%i\t%i\t%" PRIu64 "\t%" PRIu64, kmer_str, gap_pos, gap_len, background_count, signal_count);
            if (is_localmax) printf ("\tLocalmax");
            if (!only_local_max || is_localmax == 1) printf("\n");
            }
        }
    free(kmer_str);
}

int main(int argc, char* argv[]) {
    
    clock_t start = clock();
    
    if (argc != 4 && argc != 5 && argc != 6) {
        fprintf(stderr, "Usage: %s <background_file> <signal_file> <base_kmer_length> [count threshold for printing] [length_diff_cutoff]\n", argv[0]);
        fprintf(stderr, "  length_diff_cutoff: Optional parameter for local maxima detection (default: 0.35)\n");
        return 1;
    }

    char* filenames[NUM_FILES] = {argv[1], argv[2]};
    int base_kmer_length = atoi(argv[3]);
    double threshold = (argc == 5 || argc == 6) ? atof(argv[4]) : 10;             // Use command line value if provided
    double length_diff_cutoff = argc == 6 ? atof(argv[5]) : 0.35;  // Use command line value if provided
    
    // At the top of main, define these once:
    int shortest_kmer = base_kmer_length - 1;
    int longest_kmer = base_kmer_length + 1;
    int print_only_localmax = 1;
    
    // Print parameters
    // printf("DEBUG: Running with parameters:\n");
    // printf("DEBUG: Base kmer length: %d\n", base_kmer_length);
    // printf("DEBUG: Length difference cutoff: %f\n", length_diff_cutoff);
    // printf("DEBUG: Will analyze lengths %d, %d, and %d\n",
    //        base_kmer_length-1, base_kmer_length, base_kmer_length+1);
    
    GenerateMask ();
    
    // Allocate 5D results array
    uint64_t***** results = allocate_5d_results(base_kmer_length);
    if (!results) {
        fprintf(stderr, "ERROR: Memory allocation failed for results array\n");
        return 1;
    }

    // Process each file
    for (int file_idx = 0; file_idx < NUM_FILES; file_idx++) {
        LineBuffer buffer;
        buffer.results = results;
        buffer.bitstream = (uint64_t*)calloc(MAX_SEQ_LEN * 2, sizeof(uint64_t));
        
        if (!buffer.bitstream) {
            fprintf(stderr, "ERROR: Memory allocation failed for bitstream\n");
            free_5d_results(results, base_kmer_length);  // Fixed call
            return 1;
        }

        init_buffer(&buffer, base_kmer_length, file_idx, filenames[file_idx]);

        pthread_t producer_thread, consumer_thread;
        pthread_create(&producer_thread, NULL, producer, &buffer);
        pthread_create(&consumer_thread, NULL, consumer, &buffer);

        pthread_join(producer_thread, NULL);
        pthread_join(consumer_thread, NULL);

        destroy_buffer(&buffer);
        free(buffer.bitstream);
    }


    int file = 0;
    
    // PRINT KMERS
    // for (int file = 0; file < NUM_FILES; file++) {
    //    printf("\nFile %d results:\n", file);
        printf("kmer\tgap position\tgap length\tbackground\tsignal\tLocalmax\n");
        
        // First print ungapped kmers above threshold
        print_kmers_single_config(results, file, base_kmer_length, shortest_kmer, longest_kmer,
                                0, 0, threshold, print_only_localmax, length_diff_cutoff);
        
        // Then print gapped local maxima
        for (int klen = shortest_kmer; klen <= longest_kmer; klen++) {
            if (klen <= 0) continue;
            
            // Calculate center positions
            if (klen % 2 == 0) {
                // Even length: one center position
                int gap_pos = klen/2;
                for (int gap_len = 1; gap_len <= 10; gap_len++) {
                    print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                            gap_pos, gap_len, threshold, print_only_localmax, length_diff_cutoff);
                }
            } else {
                // Odd length: two center positions
                for (int gap_pos = klen/2; gap_pos <= klen/2 + 1; gap_pos++) {
                    for (int gap_len = 1; gap_len <= 10; gap_len++) {
                        print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                                gap_pos, gap_len, threshold, 1, length_diff_cutoff);
                    }
                }
            }
        }
    // }

    
    // Cleanup
    free_5d_results(results, base_kmer_length);
    
    clock_t end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total execution time: %.2f seconds\n", cpu_time_used);
    
    return 0;
}
