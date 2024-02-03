      int     FUNCTION IPARAM2STAGE( ISPEC, NAME, OPTS, NI, NBI, IBI, NXI );
#if defined(_OPENMP)
      use omp_lib
#endif
      IMPLICIT NONE

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String   *( * )    NAME, OPTS;
      int                ISPEC, NI, NBI, IBI, NXI;

*  ================================================================
      // ..
      // .. Local Scalars ..
      int                I, IC, IZ, KD, IB, LHOUS, LWORK, NTHREADS, FACTOPTNB, QROPTNB, LQOPTNB;
      bool               RPREC, CPREC;
      String             PREC*1, ALGO*3, STAG*5, SUBNAM*12, VECT*1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CHAR, ICHAR, MAX
      // ..
      // .. External Functions ..
      int                ILAENV;
      bool               LSAME;
      // EXTERNAL ILAENV, LSAME
      // ..
      // .. Executable Statements ..

      // Invalid value for ISPEC

      if ( (ISPEC.LT.17) || (ISPEC.GT.21) ) {
          IPARAM2STAGE = -1
          RETURN
      }

      // Get the number of threads

      NTHREADS = 1
#if defined(_OPENMP)
!$OMP PARALLEL
      NTHREADS = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif
       // WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC

      if ( ISPEC != 19 ) {

         // Convert NAME to upper case if the first character is lower case.

         IPARAM2STAGE = -1
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         if ( IZ == 90 || IZ == 122 ) {

            // ASCII character set

            if ( IC.GE.97 && IC.LE.122 ) {
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               for (I = 2; I <= 12; I++) { // 100
                  IC = ICHAR( SUBNAM( I: I ) )
                  if (IC.GE.97 && IC.LE.122) SUBNAM( I: I ) = CHAR( IC-32 );
               } // 100
            }

         } else if ( IZ == 233 || IZ == 169 ) {

            // EBCDIC character set

            if ( ( IC.GE.129 && IC.LE.137 ) || ( IC.GE.145 && IC.LE.153 ) || ( IC.GE.162 && IC.LE.169 ) ) {
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               for (I = 2; I <= 12; I++) { // 110
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 && IC.LE.137 ) || ( IC.GE.145 && IC.LE.153 ) || ( IC.GE.162 && IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
               } // 110
            }

         } else if ( IZ == 218 || IZ == 250 ) {

            // Prime machines:  ASCII+128

            if ( IC.GE.225 && IC.LE.250 ) {
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               for (I = 2; I <= 12; I++) { // 120
                 IC = ICHAR( SUBNAM( I: I ) )
                 if (IC.GE.225 && IC.LE.250) SUBNAM( I: I ) = CHAR( IC-32 );
               } // 120
            }
         }

         PREC  = SUBNAM( 1: 1 )
         ALGO  = SUBNAM( 4: 6 )
         STAG  = SUBNAM( 8:12 )
         RPREC = PREC == 'S' || PREC == 'D'
         CPREC = PREC == 'C' || PREC == 'Z'

         // Invalid value for PRECISION

         if ( .NOT.( RPREC || CPREC ) ) {
             IPARAM2STAGE = -1
             RETURN
         }
      }
       // WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC,
      // $           '   ALGO ',ALGO,'    STAGE ',STAG


      if (( ISPEC == 17 ) || ( ISPEC == 18 )) {

      // ISPEC = 17, 18:  block size KD, IB
      // Could be also dependent from N but for now it
      // depend only on sequential or parallel

         if ( NTHREADS.GT.4 ) {
            if ( CPREC ) {
               KD = 128
               IB = 32
            } else {
               KD = 160
               IB = 40
            }
         } else if ( NTHREADS.GT.1 ) {
            if ( CPREC ) {
               KD = 64
               IB = 32
            } else {
               KD = 64
               IB = 32
            }
         } else {
            if ( CPREC ) {
               KD = 16
               IB = 16
            } else {
               KD = 32
               IB = 16
            }
         }
         if (ISPEC == 17) IPARAM2STAGE = KD;
         if (ISPEC == 18) IPARAM2STAGE = IB;

      } else if ( ISPEC == 19 ) {

      // ISPEC = 19:
      // LHOUS length of the Houselholder representation
      // matrix (V,T) of the second stage. should be >= 1.

      // Will add the VECT OPTION HERE next release
         VECT  = OPTS(1:1)
         if ( LSAME( VECT, 'N' ) ) {
            LHOUS = MAX( 1, 4*NI )
         } else {
            // This is not correct, it need to call the ALGO and the stage2
            LHOUS = MAX( 1, 4*NI ) + IBI
         }
         if ( LHOUS.GE.0 ) {
            IPARAM2STAGE = LHOUS
         } else {
            IPARAM2STAGE = -1
         }

      } else if ( ISPEC == 20 ) {

      // ISPEC = 20: (21 for future use)
      // LWORK length of the workspace for
      // either or both stages for TRD and BRD. should be >= 1.
      // TRD:
      // TRD_stage 1: = LT + LW + LS1 + LS2
                   // = LDT*KD + N*KD + N*MAX(KD,FACTOPTNB) + LDS2*KD
                     // where LDT=LDS2=KD
                   // = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
      // TRD_stage 2: = (2NB+1)*N + KD*NTHREADS
      // TRD_both   : = max(stage1,stage2) + AB ( AB=(KD+1)*N )
                   // = N*KD + N*max(KD+1,FACTOPTNB)
                     // + max(2*KD*KD, KD*NTHREADS)
                     // + (KD+1)*N
         LWORK        = -1
         SUBNAM(1:1)  = PREC
         SUBNAM(2:6)  = 'GEQRF'
         QROPTNB      = ILAENV( 1, SUBNAM, ' ', NI, NBI, -1, -1 )
         SUBNAM(2:6)  = 'GELQF'
         LQOPTNB      = ILAENV( 1, SUBNAM, ' ', NBI, NI, -1, -1 )
         // Could be QR or LQ for TRD and the max for BRD
         FACTOPTNB    = MAX(QROPTNB, LQOPTNB)
         if ( ALGO == 'TRD' ) {
            if ( STAG == '2STAG' ) {
               LWORK = NI*NBI + NI*MAX(NBI+1,FACTOPTNB) + MAX(2*NBI*NBI, NBI*NTHREADS) + (NBI+1)*NI
            } else if ( (STAG == 'HE2HB') || (STAG == 'SY2SB') ) {
               LWORK = NI*NBI + NI*MAX(NBI,FACTOPTNB) + 2*NBI*NBI
            } else if ( (STAG == 'HB2ST') || (STAG == 'SB2ST') ) {
               LWORK = (2*NBI+1)*NI + NBI*NTHREADS
            }
         } else if ( ALGO == 'BRD' ) {
            if ( STAG == '2STAG' ) {
               LWORK = 2*NI*NBI + NI*MAX(NBI+1,FACTOPTNB) + MAX(2*NBI*NBI, NBI*NTHREADS) + (NBI+1)*NI
            } else if ( STAG == 'GE2GB' ) {
               LWORK = NI*NBI + NI*MAX(NBI,FACTOPTNB) + 2*NBI*NBI
            } else if ( STAG == 'GB2BD' ) {
               LWORK = (3*NBI+1)*NI + NBI*NTHREADS
            }
         }
         LWORK = MAX ( 1, LWORK )

         if ( LWORK.GT.0 ) {
            IPARAM2STAGE = LWORK
         } else {
            IPARAM2STAGE = -1
         }

      } else if ( ISPEC == 21 ) {

      // ISPEC = 21 for future use
         IPARAM2STAGE = NXI
      }

      // ==== End of IPARAM2STAGE ====

      }
