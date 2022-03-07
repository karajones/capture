# Preparing libraries for sequencing
This is an overview of the lab work needed to complete the entire capture process from raw DNA to finished capture library ready for sequencing. This is intended to supplement the manufacturers’ protocols with advice on what I did, especially if it pertains to making changes specific to working with salamanders. 

## Required reagents
Reagents/consumables needed to shear, library prep, and capture. I’ve listed the specific items I used but there are multiple manufacturers that make similar items that would probably work just as well. This list doesn’t include standard lab items, just the specialty stuff that’s probably not lying around already.

#### 1. Arbor BioSciences MyBaits kit
	
  - For salamanders, I’d recommend 4 individuals per reaction (i.e., a 16 Rx kit would process up to 64 individuals), but larger pools might work
	
#### 2. Multiplex Oligos for Illumina

*NEBNext® Multiplex Oligos for Illumina® (96 Unique Dual Index Primer Pairs) NEB #E6440*

#### 3. Library Prep kit for Illumina

*NEBNext® UltraTM II DNA Library Prep Kit for Illumina® NEB #E7645*
  - Whichever kit you use should include an amplification step
  - Don’t use NEXTERA unless you want to do extra work to prep for the capture reactions (see the MyBaits manual for more info)

#### 4. Covaris consumables (for Covaris M220)

*microTUBE-50 AFA Fiber Screw-Cap (PN #520166)*

#### 5. Library amplification primer mix (P5/P7 library primers)

*xGen™ Library Amplification Primer Mix #1077676*
- Any primer mix that targets the P5/P7 regions of Illumina adapters will work

#### 6. Taq Polymerase Ready Mix (for final capture amplification step)

*KAPA HiFi HotStart ReadyMix (Roche Sequencing #KK2602)*
- I tried another similar Taq ready mix and it didn’t work so I’d recommend sticking with the KAPA HiFi

#### 7. SPRI beads (for cleanup and size selection)

*AmpureXP beads (Beckman Coulter #A63881)*

## Lab procedure overview
Here’s the entire lab procedure from start to finish, with my notes. See the individual protocols provided by the manufacturers for all the necessary details.

### Step I: Shear DNA to 250-300 bp
- For salamanders, start with 2000-3000 ng DNA at this step
- I used a Covaris M220 but any *mechanical* shearing method with work
	- Do NOT use an enzyme-based shearing method!

<img src="images/2100 expert_High Sensitivity DNA Assay_DE13804202_2019-03-11_17-09-37_EGRAM_Sample2.bmp">

>Bioanalyzer trace for an individual after shearing, ready to go into library prep.


### Step II: Library Prep
Library prep adds indexes and size selects the DNA to a tighter range (which improves capture/sequencing efficiency). All reagents included in NEBNext UltraTM II DNA Library Prep Kit for Illumina except where noted with *required* tag.
	
#### 1. End prep

  - I use the full volume from shearing (which is a little more than the recommended 50 ul)
	
#### 2. Adapter ligation

*Required*: Adaptor for Illumina and USER enzyme (both included with NEBNext Multiplex Oligos for Illumina)
  - Do not dilute the adapters
  
#### 3. Size selection

*Required*: SPRI beads
  - Use size selection recommendations in manual for 250 bp insert size
  
#### 4. PCR Amplification

*Required*: NEBNext Multiplex Oligos for Illumina
  - If using the NEBNext oligos, the forward and reverse primer are already combined
  - For salamanders, use 8 cycles of amplification
    
#### 5. Bead cleanup
*Required*: SPRI beads (same as size selection step)
  - You can pool individuals that are going to be in the same capture reaction before this bead cleanup step and clean them all together
  - Make sure to adjust the bead ratio to account for the larger volume of pooled individuals
  - After bead cleanup, use a vacufuge (or similar) to concentrate the library down to 7 μl before continuing with bait capture

<img src="images/2100 expert_High Sensitivity DNA Assay_DE13804202_2019-03-15_15-03-56_EGRAM_Sample1.bmp">

>Bioanalyzer trace after library prep with size selection for 250 bp. Keep in mind that a 250 bp insert size translates to a 370 bp length fragment (since the adapters and primers add 120 bp to each fragment).

## Step III: Capture using MyBaits (standard protocol)	

All reagents included in Arbor BioSciences MyBaits kit except where noted with *required* tag.

#### 1. Hybridization
For salamanders, I’d recommend the longer 24 hour hybridization time

#### 2. Bind and wash cleanup
No changes from standard protocol

#### 3. Capture library amplification

*Required*: 2X KAPA HiFi HotStart Ready Mix *and* Library amplification primer mix
  - Recommend 8 cycles of amplification (but see notes below bioanalyzer traces)

#### 4. Final bead cleanup

*Required:* SPRI beads
  - Libraries can be further pooled before bead cleanup
  - Recommend a 2:1 or 1.8:1 ratio of beads:reaction volume

<img src="images/2100 expert_High Sensitivity DNA Assay_DE13804202_2019-05-01_10-52-02_EGRAM_Sample1.bmp">
<img src="images/2100 expert_High Sensitivity DNA Assay_DE13804202_2019-05-01_10-52-02_EGRAM_Sample2.bmp">

>Bioanalyzer traces after capture and amplification. Top: finished capture library diluted before bioanalyzing. Bottom: same library undiluted.

### Recommendations on final amplification
This library was overamplified (14 cycles of amplification, over 600 ng/μl final concentration). About 45% of the mapped reads had to be removed as duplicates. :persevere: Moral of the story: don’t overamplify! If you start out with the same amount of DNA for each individual then test one capture library (4 pooled individuals) at 8 cycles and quantify it. Usually sequencing centers require about 2 ng/μl concentration, so aim for that or a little above. Adjust the number of cycles based on your test amplification. This will also reduce the number of unwanted fragments in the final library (i.e., the hump of fragments over 500 bp you can see on the second bioanalyzer trace above).
