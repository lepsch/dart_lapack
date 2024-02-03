      SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ;
      int                IFST, ILST, INFO, LDQ, LDT, N;
      // ..
      // .. Array Arguments ..
      double             Q( LDQ, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
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
      // EXTERNAL DLAEXC, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input arguments.

      INFO = 0;
      WANTQ = LSAME( COMPQ, 'V' );
      if ( !WANTQ && !LSAME( COMPQ, 'N' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDT < MAX( 1, N ) ) {
         INFO = -4;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < MAX( 1, N ) ) ) {
         INFO = -6;
      } else if (( IFST < 1 || IFST > N ) && ( N > 0 )) {
         INFO = -7;
      } else if (( ILST < 1 || ILST > N ) && ( N > 0 )) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('DTREXC', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 1) RETURN;

      // Determine the first row of specified block
      // and find out it is 1 by 1 or 2 by 2.

      if ( IFST > 1 ) {
         if( T( IFST, IFST-1 ) != ZERO ) IFST = IFST - 1;
      }
      NBF = 1;
      if ( IFST < N ) {
         if( T( IFST+1, IFST ) != ZERO ) NBF = 2;
      }

      // Determine the first row of the final block
      // and find out it is 1 by 1 or 2 by 2.

      if ( ILST > 1 ) {
         if( T( ILST, ILST-1 ) != ZERO ) ILST = ILST - 1;
      }
      NBL = 1;
      if ( ILST < N ) {
         if( T( ILST+1, ILST ) != ZERO ) NBL = 2;
      }

      if (IFST == ILST) RETURN;

      if ( IFST < ILST ) {

         // Update ILST

         if (NBF == 2 && NBL == 1) ILST = ILST - 1;
         IF( NBF == 1 && NBL == 2 ) ILST = ILST + 1;

         HERE = IFST;

         } // 10

         // Swap block with next one below

         if ( NBF == 1 || NBF == 2 ) {

            // Current block either 1 by 1 or 2 by 2

            NBNEXT = 1;
            if ( HERE+NBF+1 <= N ) {
               if( T( HERE+NBF+1, HERE+NBF ) != ZERO ) NBNEXT = 2;
            }
            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, WORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE;
               return;
            }
            HERE = HERE + NBNEXT;

            // Test if 2 by 2 block breaks into two 1 by 1 blocks

            if ( NBF == 2 ) {
               if( T( HERE+1, HERE ) == ZERO ) NBF = 3;
            }

         } else {

            // Current block consists of two 1 by 1 blocks each of which
            // must be swapped individually

            NBNEXT = 1;
            if ( HERE+3 <= N ) {
               if( T( HERE+3, HERE+2 ) != ZERO ) NBNEXT = 2;
            }
            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT, WORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE;
               return;
            }
            if ( NBNEXT == 1 ) {

               // Swap two 1 by 1 blocks, no problems possible

               dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, WORK, INFO );
               HERE = HERE + 1;
            } else {

               // Recompute NBNEXT in case 2 by 2 split

               if( T( HERE+2, HERE+1 ) == ZERO ) NBNEXT = 1;
               if ( NBNEXT == 2 ) {

                  // 2 by 2 Block did not split

                  dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, WORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE;
                     return;
                  }
                  HERE = HERE + 2;
               } else {

                  // 2 by 2 Block did split

                  dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, WORK, INFO );
                  dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1, WORK, INFO );
                  HERE = HERE + 2;
               }
            }
         }
         if (HERE < ILST) GO TO 10;

      } else {

         HERE = IFST;
         } // 20

         // Swap block with next one above

         if ( NBF == 1 || NBF == 2 ) {

            // Current block either 1 by 1 or 2 by 2

            NBNEXT = 1;
            if ( HERE >= 3 ) {
               if( T( HERE-1, HERE-2 ) != ZERO ) NBNEXT = 2;
            }
            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, NBF, WORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE;
               return;
            }
            HERE = HERE - NBNEXT;

            // Test if 2 by 2 block breaks into two 1 by 1 blocks

            if ( NBF == 2 ) {
               if( T( HERE+1, HERE ) == ZERO ) NBF = 3;
            }

         } else {

            // Current block consists of two 1 by 1 blocks each of which
            // must be swapped individually

            NBNEXT = 1;
            if ( HERE >= 3 ) {
               if( T( HERE-1, HERE-2 ) != ZERO ) NBNEXT = 2;
            }
            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, 1, WORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE;
               return;
            }
            if ( NBNEXT == 1 ) {

               // Swap two 1 by 1 blocks, no problems possible

               dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, WORK, INFO );
               HERE = HERE - 1;
            } else {

               // Recompute NBNEXT in case 2 by 2 split

               if( T( HERE, HERE-1 ) == ZERO ) NBNEXT = 1;
               if ( NBNEXT == 2 ) {

                  // 2 by 2 Block did not split

                  dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1, WORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE;
                     return;
                  }
                  HERE = HERE - 2;
               } else {

                  // 2 by 2 Block did split

                  dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, WORK, INFO );
                  dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1, WORK, INFO );
                  HERE = HERE - 2;
               }
            }
         }
         if (HERE > ILST) GO TO 20;
      }
      ILST = HERE;

      return;

      // End of DTREXC

      }
