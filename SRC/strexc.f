      SUBROUTINE STREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ;
      int                IFST, ILST, INFO, LDQ, LDT, N;
      // ..
      // .. Array Arguments ..
      REAL               Q( LDQ, * ), T( LDT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               WANTQ;
      int                HERE, NBF, NBL, NBNEXT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAEXC, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input arguments.

      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      if ( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) {
         INFO = -6
      } else if (( IFST.LT.1 .OR. IFST.GT.N ).AND.( N.GT.0 )) {
         INFO = -7
      } else if (( ILST.LT.1 .OR. ILST.GT.N ).AND.( N.GT.0 )) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('STREXC', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.LE.1 ) RETURN

      // Determine the first row of specified block
      // and find out it is 1 by 1 or 2 by 2.

      if ( IFST.GT.1 ) {
         IF( T( IFST, IFST-1 ).NE.ZERO ) IFST = IFST - 1
      }
      NBF = 1
      if ( IFST.LT.N ) {
         IF( T( IFST+1, IFST ).NE.ZERO ) NBF = 2
      }

      // Determine the first row of the final block
      // and find out it is 1 by 1 or 2 by 2.

      if ( ILST.GT.1 ) {
         IF( T( ILST, ILST-1 ).NE.ZERO ) ILST = ILST - 1
      }
      NBL = 1
      if ( ILST.LT.N ) {
         IF( T( ILST+1, ILST ).NE.ZERO ) NBL = 2
      }

      IF( IFST.EQ.ILST ) RETURN

      if ( IFST.LT.ILST ) {

         // Update ILST

         IF( NBF.EQ.2 .AND. NBL.EQ.1 ) ILST = ILST - 1          IF( NBF.EQ.1 .AND. NBL.EQ.2 ) ILST = ILST + 1

         HERE = IFST

         } // 10

         // Swap block with next one below

         if ( NBF.EQ.1 .OR. NBF.EQ.2 ) {

            // Current block either 1 by 1 or 2 by 2

            NBNEXT = 1
            if ( HERE+NBF+1.LE.N ) {
               IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO ) NBNEXT = 2
            }
            slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, WORK, INFO );
            if ( INFO.NE.0 ) {
               ILST = HERE
               RETURN
            }
            HERE = HERE + NBNEXT

            // Test if 2 by 2 block breaks into two 1 by 1 blocks

            if ( NBF.EQ.2 ) {
               IF( T( HERE+1, HERE ).EQ.ZERO ) NBF = 3
            }

         } else {

            // Current block consists of two 1 by 1 blocks each of which
            // must be swapped individually

            NBNEXT = 1
            if ( HERE+3.LE.N ) {
               IF( T( HERE+3, HERE+2 ).NE.ZERO ) NBNEXT = 2
            }
            slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT, WORK, INFO );
            if ( INFO.NE.0 ) {
               ILST = HERE
               RETURN
            }
            if ( NBNEXT.EQ.1 ) {

               // Swap two 1 by 1 blocks, no problems possible

               slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, WORK, INFO );
               HERE = HERE + 1
            } else {

               // Recompute NBNEXT in case 2 by 2 split

               IF( T( HERE+2, HERE+1 ).EQ.ZERO ) NBNEXT = 1
               if ( NBNEXT.EQ.2 ) {

                  // 2 by 2 Block did not split

                  slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, WORK, INFO );
                  if ( INFO.NE.0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE + 2
               } else {

                  // 2 by 2 Block did split

                  slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, WORK, INFO )                   CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1, WORK, INFO );
                  HERE = HERE + 2
               }
            }
         }
         IF( HERE.LT.ILST ) GO TO 10

      } else {

         HERE = IFST
         } // 20

         // Swap block with next one above

         if ( NBF.EQ.1 .OR. NBF.EQ.2 ) {

            // Current block either 1 by 1 or 2 by 2

            NBNEXT = 1
            if ( HERE.GE.3 ) {
               IF( T( HERE-1, HERE-2 ).NE.ZERO ) NBNEXT = 2
            }
            slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, NBF, WORK, INFO );
            if ( INFO.NE.0 ) {
               ILST = HERE
               RETURN
            }
            HERE = HERE - NBNEXT

            // Test if 2 by 2 block breaks into two 1 by 1 blocks

            if ( NBF.EQ.2 ) {
               IF( T( HERE+1, HERE ).EQ.ZERO ) NBF = 3
            }

         } else {

            // Current block consists of two 1 by 1 blocks each of which
            // must be swapped individually

            NBNEXT = 1
            if ( HERE.GE.3 ) {
               IF( T( HERE-1, HERE-2 ).NE.ZERO ) NBNEXT = 2
            }
            slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, 1, WORK, INFO );
            if ( INFO.NE.0 ) {
               ILST = HERE
               RETURN
            }
            if ( NBNEXT.EQ.1 ) {

               // Swap two 1 by 1 blocks, no problems possible

               slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, WORK, INFO );
               HERE = HERE - 1
            } else {

               // Recompute NBNEXT in case 2 by 2 split

               IF( T( HERE, HERE-1 ).EQ.ZERO ) NBNEXT = 1
               if ( NBNEXT.EQ.2 ) {

                  // 2 by 2 Block did not split

                  slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1, WORK, INFO );
                  if ( INFO.NE.0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE - 2
               } else {

                  // 2 by 2 Block did split

                  slaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, WORK, INFO )                   CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1, WORK, INFO );
                  HERE = HERE - 2
               }
            }
         }
         IF( HERE.GT.ILST ) GO TO 20
      }
      ILST = HERE

      RETURN

      // End of STREXC

      }
