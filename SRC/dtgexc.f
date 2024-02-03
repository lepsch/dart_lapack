      SUBROUTINE DTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IFST, ILST, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                HERE, LWMIN, NBF, NBL, NBNEXT;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTGEX2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Decode and test input arguments.

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDQ.LT.1 || WANTQ && ( LDQ.LT.MAX( 1, N ) ) ) {
         INFO = -9
      } else if ( LDZ.LT.1 || WANTZ && ( LDZ.LT.MAX( 1, N ) ) ) {
         INFO = -11
      } else if ( IFST.LT.1 || IFST.GT.N ) {
         INFO = -12
      } else if ( ILST.LT.1 || ILST.GT.N ) {
         INFO = -13
      }

      if ( INFO == 0 ) {
         if ( N.LE.1 ) {
            LWMIN = 1
         } else {
            LWMIN = 4*N + 16
         }
         WORK(1) = LWMIN

         if (LWORK.LT.LWMIN && .NOT.LQUERY) {
            INFO = -15
         }
      }

      if ( INFO != 0 ) {
         xerbla('DTGEXC', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.LE.1) RETURN;

      // Determine the first row of the specified block and find out
      // if it is 1-by-1 or 2-by-2.

      if ( IFST.GT.1 ) {
         IF( A( IFST, IFST-1 ) != ZERO ) IFST = IFST - 1
      }
      NBF = 1
      if ( IFST.LT.N ) {
         IF( A( IFST+1, IFST ) != ZERO ) NBF = 2
      }

      // Determine the first row of the final block
      // and find out if it is 1-by-1 or 2-by-2.

      if ( ILST.GT.1 ) {
         IF( A( ILST, ILST-1 ) != ZERO ) ILST = ILST - 1
      }
      NBL = 1
      if ( ILST.LT.N ) {
         IF( A( ILST+1, ILST ) != ZERO ) NBL = 2
      }
      if (IFST == ILST) RETURN;

      if ( IFST.LT.ILST ) {

         // Update ILST.

         if (NBF == 2 && NBL == 1) ILST = ILST - 1          IF( NBF == 1 && NBL == 2 ) ILST = ILST + 1;

         HERE = IFST

         } // 10

         // Swap with next one below.

         if ( NBF == 1 || NBF == 2 ) {

            // Current block either 1-by-1 or 2-by-2.

            NBNEXT = 1
            if ( HERE+NBF+1.LE.N ) {
               IF( A( HERE+NBF+1, HERE+NBF ) != ZERO ) NBNEXT = 2
            }
            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, NBF, NBNEXT, WORK, LWORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE
               RETURN
            }
            HERE = HERE + NBNEXT

            // Test if 2-by-2 block breaks into two 1-by-1 blocks.

            if ( NBF == 2 ) {
               IF( A( HERE+1, HERE ) == ZERO ) NBF = 3
            }

         } else {

            // Current block consists of two 1-by-1 blocks, each of which
            // must be swapped individually.

            NBNEXT = 1
            if ( HERE+3.LE.N ) {
               IF( A( HERE+3, HERE+2 ) != ZERO ) NBNEXT = 2
            }
            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE+1, 1, NBNEXT, WORK, LWORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE
               RETURN
            }
            if ( NBNEXT == 1 ) {

               // Swap two 1-by-1 blocks.

               dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO );
               if ( INFO != 0 ) {
                  ILST = HERE
                  RETURN
               }
               HERE = HERE + 1

            } else {

               // Recompute NBNEXT in case of 2-by-2 split.

               IF( A( HERE+2, HERE+1 ) == ZERO ) NBNEXT = 1
               if ( NBNEXT == 2 ) {

                  // 2-by-2 block did not split.

                  dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, NBNEXT, WORK, LWORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE + 2
               } else {

                  // 2-by-2 block did split.

                  dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE + 1
                  dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE + 1
               }

            }
         }
         if (HERE.LT.ILST) GO TO 10;
      } else {
         HERE = IFST

         } // 20

         // Swap with next one below.

         if ( NBF == 1 || NBF == 2 ) {

            // Current block either 1-by-1 or 2-by-2.

            NBNEXT = 1
            if ( HERE.GE.3 ) {
               IF( A( HERE-1, HERE-2 ) != ZERO ) NBNEXT = 2
            }
            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE-NBNEXT, NBNEXT, NBF, WORK, LWORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE
               RETURN
            }
            HERE = HERE - NBNEXT

            // Test if 2-by-2 block breaks into two 1-by-1 blocks.

            if ( NBF == 2 ) {
               IF( A( HERE+1, HERE ) == ZERO ) NBF = 3
            }

         } else {

            // Current block consists of two 1-by-1 blocks, each of which
            // must be swapped individually.

            NBNEXT = 1
            if ( HERE.GE.3 ) {
               IF( A( HERE-1, HERE-2 ) != ZERO ) NBNEXT = 2
            }
            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE-NBNEXT, NBNEXT, 1, WORK, LWORK, INFO );
            if ( INFO != 0 ) {
               ILST = HERE
               RETURN
            }
            if ( NBNEXT == 1 ) {

               // Swap two 1-by-1 blocks.

               dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, NBNEXT, 1, WORK, LWORK, INFO );
               if ( INFO != 0 ) {
                  ILST = HERE
                  RETURN
               }
               HERE = HERE - 1
            } else {

              // Recompute NBNEXT in case of 2-by-2 split.

               IF( A( HERE, HERE-1 ) == ZERO ) NBNEXT = 1
               if ( NBNEXT == 2 ) {

                  // 2-by-2 block did not split.

                  dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE-1, 2, 1, WORK, LWORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE - 2
               } else {

                  // 2-by-2 block did split.

                  dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE - 1
                  dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO );
                  if ( INFO != 0 ) {
                     ILST = HERE
                     RETURN
                  }
                  HERE = HERE - 1
               }
            }
         }
         if (HERE.GT.ILST) GO TO 20;
      }
      ILST = HERE
      WORK( 1 ) = LWMIN
      RETURN

      // End of DTGEXC

      }
