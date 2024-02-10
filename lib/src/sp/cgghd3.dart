      void cgghd3(COMPQ, COMPZ, N, ILO, IHI, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> Q, final int LDQ, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--


      String             COMPQ, COMPZ;
      int                IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK;
      Complex            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), WORK( * );
      // ..

      Complex            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      bool               BLK22, INITQ, INITZ, LQUERY, WANTQ, WANTZ;
      String             COMPQ2, COMPZ2;
      int                COLA, I, IERR, J, J0, JCOL, JJ, JROW, K, KACC22, LEN, LWKOPT, N2NB, NB, NBLST, NBMIN, NH, NNB, NX, PPW, PPWO, PW, TOP, TOPQ;
      double               C;
      Complex            C1, C2, CTEMP, S, S1, S2, TEMP, TEMP1, TEMP2, TEMP3;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGGHRD, CLARTG, CLASET, CUNM22, CROT, CGEMM, CGEMV, CTRMV, CLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, CMPLX, CONJG, MAX

      // Decode and test the input parameters.

      INFO = 0;
      NB = ilaenv( 1, 'CGGHD3', ' ', N, ILO, IHI, -1 );
      NH = IHI - ILO + 1;
      if ( NH <= 1 ) {
         LWKOPT = 1;
      } else {
         LWKOPT = 6*N*NB;
      }
      WORK[1] = SROUNDUP_LWORK( LWKOPT );
      INITQ = lsame( COMPQ, 'I' );
      WANTQ = INITQ || lsame( COMPQ, 'V' );
      INITZ = lsame( COMPZ, 'I' );
      WANTZ = INITZ || lsame( COMPZ, 'V' );
      LQUERY = ( LWORK == -1 );

      if ( !lsame( COMPQ, 'N' ) && !WANTQ ) {
         INFO = -1;
      } else if ( !lsame( COMPZ, 'N' ) && !WANTZ ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 ) {
         INFO = -4;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( ( WANTQ && LDQ < N ) || LDQ < 1 ) {
         INFO = -11;
      } else if ( ( WANTZ && LDZ < N ) || LDZ < 1 ) {
         INFO = -13;
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -15;
      }
      if ( INFO != 0 ) {
         xerbla('CGGHD3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Initialize Q and Z if desired.

      if (INITQ) claset( 'All', N, N, CZERO, CONE, Q, LDQ );
      IF( INITZ ) claset( 'All', N, N, CZERO, CONE, Z, LDZ );

      // Zero out lower triangle of B.

      if (N > 1) claset( 'Lower', N-1, N-1, CZERO, CZERO, B(2, 1), LDB );

      // Quick return if possible

      if ( NH <= 1 ) {
         WORK[1] = CONE;
         return;
      }

      // Determine the blocksize.

      NBMIN = ilaenv( 2, 'CGGHD3', ' ', N, ILO, IHI, -1 );
      if ( NB > 1 && NB < NH ) {

         // Determine when to use unblocked instead of blocked code.

         NX = max( NB, ilaenv( 3, 'CGGHD3', ' ', N, ILO, IHI, -1 ) );
         if ( NX < NH ) {

            // Determine if workspace is large enough for blocked code.

            if ( LWORK < LWKOPT ) {

               // Not enough workspace to use optimal NB:  determine the
               // minimum value of NB, and reduce NB or force use of
               // unblocked code.

               NBMIN = max( 2, ilaenv( 2, 'CGGHD3', ' ', N, ILO, IHI, -1 ) );
               if ( LWORK >= 6*N*NBMIN ) {
                  NB = LWORK / ( 6*N );
               } else {
                  NB = 1;
               }
            }
         }
      }

      if ( NB < NBMIN || NB >= NH ) {

         // Use unblocked code below

         JCOL = ILO;

      } else {

         // Use blocked code

         KACC22 = ilaenv( 16, 'CGGHD3', ' ', N, ILO, IHI, -1 );
         BLK22 = KACC22 == 2;
         for (JCOL = ILO; NB < 0 ? JCOL >= IHI-2 : JCOL <= IHI-2; JCOL += NB) {
            NNB = min( NB, IHI-JCOL-1 );

            // Initialize small unitary factors that will hold the
            // accumulated Givens rotations in workspace.
            // N2NB   denotes the number of 2*NNB-by-2*NNB factors
            // NBLST  denotes the (possibly smaller) order of the last
            //        factor.

            N2NB = ( IHI-JCOL-1 ) / NNB - 1;
            NBLST = IHI - JCOL - N2NB*NNB;
            claset('All', NBLST, NBLST, CZERO, CONE, WORK, NBLST );
            PW = NBLST * NBLST + 1;
            for (I = 1; I <= N2NB; I++) {
               claset('All', 2*NNB, 2*NNB, CZERO, CONE, WORK( PW ), 2*NNB );
               PW = PW + 4*NNB*NNB;
            }

            // Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form.

            for (J = JCOL; J <= JCOL+NNB-1; J++) {

               // Reduce Jth column of A. Store cosines and sines in Jth
               // column of A and B, respectively.

               for (I = IHI; I >= J+2; I--) {
                  TEMP = A( I-1, J );
                  clartg(TEMP, A( I, J ), C, S, A( I-1, J ) );
                  A[I][J] = CMPLX( C );
                  B[I][J] = S;
               }

               // Accumulate Givens rotations into workspace array.

               PPW  = ( NBLST + 1 )*( NBLST - 2 ) - J + JCOL + 1;
               LEN  = 2 + J - JCOL;
               JROW = J + N2NB*NNB + 2;
               for (I = IHI; I >= JROW; I--) {
                  CTEMP = A( I, J );
                  S = B( I, J );
                  for (JJ = PPW; JJ <= PPW+LEN-1; JJ++) {
                     TEMP = WORK( JJ + NBLST );
                     WORK[JJ + NBLST] = CTEMP*TEMP - S*WORK( JJ );
                     WORK[JJ] = CONJG( S )*TEMP + CTEMP*WORK( JJ );
                  }
                  LEN = LEN + 1;
                  PPW = PPW - NBLST - 1;
               }

               PPWO = NBLST*NBLST + ( NNB+J-JCOL-1 )*2*NNB + NNB;
               J0 = JROW - NNB;
               for (JROW = J0; -NNB < 0 ? JROW >= J+2 : JROW <= J+2; JROW += -NNB) {
                  PPW = PPWO;
                  LEN  = 2 + J - JCOL;
                  for (I = JROW+NNB-1; I >= JROW; I--) {
                     CTEMP = A( I, J );
                     S = B( I, J );
                     for (JJ = PPW; JJ <= PPW+LEN-1; JJ++) {
                        TEMP = WORK( JJ + 2*NNB );
                        WORK[JJ + 2*NNB] = CTEMP*TEMP - S*WORK( JJ );
                        WORK[JJ] = CONJG( S )*TEMP + CTEMP*WORK( JJ );
                     }
                     LEN = LEN + 1;
                     PPW = PPW - 2*NNB - 1;
                  }
                  PPWO = PPWO + 4*NNB*NNB;
               }

               // TOP denotes the number of top rows in A and B that will
               // not be updated during the next steps.

               if ( JCOL <= 2 ) {
                  TOP = 0;
               } else {
                  TOP = JCOL;
               }

               // Propagate transformations through B and replace stored
               // left sines/cosines by right sines/cosines.

               for (JJ = N; JJ >= J+1; JJ--) {

                  // Update JJth column of B.

                  for (I = min( JJ+1, IHI ); I >= J+2; I--) {
                     CTEMP = A( I, J );
                     S = B( I, J );
                     TEMP = B( I, JJ );
                     B[I][JJ] = CTEMP*TEMP - CONJG( S )*B( I-1, JJ );
                     B[I-1][JJ] = S*TEMP + CTEMP*B( I-1, JJ );
                  }

                  // Annihilate B( JJ+1, JJ ).

                  if ( JJ < IHI ) {
                     TEMP = B( JJ+1, JJ+1 );
                     clartg(TEMP, B( JJ+1, JJ ), C, S, B( JJ+1, JJ+1 ) );
                     B[JJ+1][JJ] = CZERO;
                     crot(JJ-TOP, B( TOP+1, JJ+1 ), 1, B( TOP+1, JJ ), 1, C, S );
                     A[JJ+1][J] = CMPLX( C );
                     B[JJ+1][J] = -CONJG( S );
                  }
               }

               // Update A by transformations from right.

               JJ = (IHI-J-1 % 3);
               for (I = IHI-J-3; I >= JJ+1; I -= 3) {
                  CTEMP = A( J+1+I, J );
                  S = -B( J+1+I, J );
                  C1 = A( J+2+I, J );
                  S1 = -B( J+2+I, J );
                  C2 = A( J+3+I, J );
                  S2 = -B( J+3+I, J );

                  for (K = TOP+1; K <= IHI; K++) {
                     TEMP = A( K, J+I  );
                     TEMP1 = A( K, J+I+1 );
                     TEMP2 = A( K, J+I+2 );
                     TEMP3 = A( K, J+I+3 );
                     A[K][J+I+3] = C2*TEMP3 + CONJG( S2 )*TEMP2;
                     TEMP2 = -S2*TEMP3 + C2*TEMP2;
                     A[K][J+I+2] = C1*TEMP2 + CONJG( S1 )*TEMP1;
                     TEMP1 = -S1*TEMP2 + C1*TEMP1;
                     A[K][J+I+1] = CTEMP*TEMP1 + CONJG( S )*TEMP;
                     A[K][J+I] = -S*TEMP1 + CTEMP*TEMP;
                  }
               }

               if ( JJ > 0 ) {
                  for (I = JJ; I >= 1; I--) {
                     C = double( A( J+1+I, J ) );
                     crot(IHI-TOP, A( TOP+1, J+I+1 ), 1, A( TOP+1, J+I ), 1, C, -CONJG( B( J+1+I, J ) ) );
                  }
               }

               // Update (J+1)th column of A by transformations from left.

               if ( J < JCOL + NNB - 1 ) {
                  LEN  = 1 + J - JCOL;

                  // Multiply with the trailing accumulated unitary
                  // matrix, which takes the form

                         // [  U11  U12  ]
                     // U = [            ],
                     //     [  U21  U22  ]

                  // where U21 is a LEN-by-LEN matrix and U12 is lower
                  // triangular.

                  JROW = IHI - NBLST + 1;
                  cgemv('Conjugate', NBLST, LEN, CONE, WORK, NBLST, A( JROW, J+1 ), 1, CZERO, WORK( PW ), 1 );
                  PPW = PW + LEN;
                  for (I = JROW; I <= JROW+NBLST-LEN-1; I++) {
                     WORK[PPW] = A( I, J+1 );
                     PPW = PPW + 1;
                  }
                  ctrmv('Lower', 'Conjugate', 'Non-unit', NBLST-LEN, WORK( LEN*NBLST + 1 ), NBLST, WORK( PW+LEN ), 1 );
                  cgemv('Conjugate', LEN, NBLST-LEN, CONE, WORK( (LEN+1)*NBLST - LEN + 1 ), NBLST, A( JROW+NBLST-LEN, J+1 ), 1, CONE, WORK( PW+LEN ), 1 );
                  PPW = PW;
                  for (I = JROW; I <= JROW+NBLST-1; I++) {
                     A[I][J+1] = WORK( PPW );
                     PPW = PPW + 1;
                  }

                  // Multiply with the other accumulated unitary
                  // matrices, which take the form

                         // [  U11  U12   0  ]
                         // [                ]
                     // U = [  U21  U22   0  ],
                     //     [                ]
                     //     [   0    0    I  ]

                  // where I denotes the (NNB-LEN)-by-(NNB-LEN) identity
                  // matrix, U21 is a LEN-by-LEN upper triangular matrix
                  // and U12 is an NNB-by-NNB lower triangular matrix.

                  PPWO = 1 + NBLST*NBLST;
                  J0 = JROW - NNB;
                  for (JROW = J0; -NNB < 0 ? JROW >= JCOL+1 : JROW <= JCOL+1; JROW += -NNB) {
                     PPW = PW + LEN;
                     for (I = JROW; I <= JROW+NNB-1; I++) {
                        WORK[PPW] = A( I, J+1 );
                        PPW = PPW + 1;
                     }
                     PPW = PW;
                     for (I = JROW+NNB; I <= JROW+NNB+LEN-1; I++) {
                        WORK[PPW] = A( I, J+1 );
                        PPW = PPW + 1;
                     }
                     ctrmv('Upper', 'Conjugate', 'Non-unit', LEN, WORK( PPWO + NNB ), 2*NNB, WORK( PW ), 1 );
                     ctrmv('Lower', 'Conjugate', 'Non-unit', NNB, WORK( PPWO + 2*LEN*NNB ), 2*NNB, WORK( PW + LEN ), 1 );
                     cgemv('Conjugate', NNB, LEN, CONE, WORK( PPWO ), 2*NNB, A( JROW, J+1 ), 1, CONE, WORK( PW ), 1 );
                     cgemv('Conjugate', LEN, NNB, CONE, WORK( PPWO + 2*LEN*NNB + NNB ), 2*NNB, A( JROW+NNB, J+1 ), 1, CONE, WORK( PW+LEN ), 1 );
                     PPW = PW;
                     for (I = JROW; I <= JROW+LEN+NNB-1; I++) {
                        A[I][J+1] = WORK( PPW );
                        PPW = PPW + 1;
                     }
                     PPWO = PPWO + 4*NNB*NNB;
                  }
               }
            }

            // Apply accumulated unitary matrices to A.

            COLA = N - JCOL - NNB + 1;
            J = IHI - NBLST + 1;
            cgemm('Conjugate', 'No Transpose', NBLST, COLA, NBLST, CONE, WORK, NBLST, A( J, JCOL+NNB ), LDA, CZERO, WORK( PW ), NBLST );
            clacpy('All', NBLST, COLA, WORK( PW ), NBLST, A( J, JCOL+NNB ), LDA );
            PPWO = NBLST*NBLST + 1;
            J0 = J - NNB;
            for (J = J0; -NNB < 0 ? J >= JCOL+1 : J <= JCOL+1; J += -NNB) {
               if ( BLK22 ) {

                  // Exploit the structure of

                         // [  U11  U12  ]
                     // U = [            ]
                     //     [  U21  U22  ],

                  // where all blocks are NNB-by-NNB, U21 is upper
                  // triangular and U12 is lower triangular.

                  cunm22('Left', 'Conjugate', 2*NNB, COLA, NNB, NNB, WORK( PPWO ), 2*NNB, A( J, JCOL+NNB ), LDA, WORK( PW ), LWORK-PW+1, IERR );
               } else {

                  // Ignore the structure of U.

                  cgemm('Conjugate', 'No Transpose', 2*NNB, COLA, 2*NNB, CONE, WORK( PPWO ), 2*NNB, A( J, JCOL+NNB ), LDA, CZERO, WORK( PW ), 2*NNB );
                  clacpy('All', 2*NNB, COLA, WORK( PW ), 2*NNB, A( J, JCOL+NNB ), LDA );
               }
               PPWO = PPWO + 4*NNB*NNB;
            }

            // Apply accumulated unitary matrices to Q.

            if ( WANTQ ) {
               J = IHI - NBLST + 1;
               if ( INITQ ) {
                  TOPQ = max( 2, J - JCOL + 1 );
                  NH  = IHI - TOPQ + 1;
               } else {
                  TOPQ = 1;
                  NH = N;
               }
               cgemm('No Transpose', 'No Transpose', NH, NBLST, NBLST, CONE, Q( TOPQ, J ), LDQ, WORK, NBLST, CZERO, WORK( PW ), NH );
               clacpy('All', NH, NBLST, WORK( PW ), NH, Q( TOPQ, J ), LDQ );
               PPWO = NBLST*NBLST + 1;
               J0 = J - NNB;
               for (J = J0; -NNB < 0 ? J >= JCOL+1 : J <= JCOL+1; J += -NNB) {
                  if ( INITQ ) {
                     TOPQ = max( 2, J - JCOL + 1 );
                     NH  = IHI - TOPQ + 1;
                  }
                  if ( BLK22 ) {

                     // Exploit the structure of U.

                     cunm22('Right', 'No Transpose', NH, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, Q( TOPQ, J ), LDQ, WORK( PW ), LWORK-PW+1, IERR );
                  } else {

                     // Ignore the structure of U.

                     cgemm('No Transpose', 'No Transpose', NH, 2*NNB, 2*NNB, CONE, Q( TOPQ, J ), LDQ, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), NH );
                     clacpy('All', NH, 2*NNB, WORK( PW ), NH, Q( TOPQ, J ), LDQ );
                  }
                  PPWO = PPWO + 4*NNB*NNB;
               }
            }

            // Accumulate right Givens rotations if required.

            if ( WANTZ || TOP > 0 ) {

               // Initialize small unitary factors that will hold the
               // accumulated Givens rotations in workspace.

               claset('All', NBLST, NBLST, CZERO, CONE, WORK, NBLST );
               PW = NBLST * NBLST + 1;
               for (I = 1; I <= N2NB; I++) {
                  claset('All', 2*NNB, 2*NNB, CZERO, CONE, WORK( PW ), 2*NNB );
                  PW = PW + 4*NNB*NNB;
               }

               // Accumulate Givens rotations into workspace array.

               for (J = JCOL; J <= JCOL+NNB-1; J++) {
                  PPW  = ( NBLST + 1 )*( NBLST - 2 ) - J + JCOL + 1;
                  LEN  = 2 + J - JCOL;
                  JROW = J + N2NB*NNB + 2;
                  for (I = IHI; I >= JROW; I--) {
                     CTEMP = A( I, J );
                     A[I][J] = CZERO;
                     S = B( I, J );
                     B[I][J] = CZERO;
                     for (JJ = PPW; JJ <= PPW+LEN-1; JJ++) {
                        TEMP = WORK( JJ + NBLST );
                        WORK[JJ + NBLST] = CTEMP*TEMP - CONJG( S )*WORK( JJ );
                        WORK[JJ] = S*TEMP + CTEMP*WORK( JJ );
                     }
                     LEN = LEN + 1;
                     PPW = PPW - NBLST - 1;
                  }

                  PPWO = NBLST*NBLST + ( NNB+J-JCOL-1 )*2*NNB + NNB;
                  J0 = JROW - NNB;
                  for (JROW = J0; -NNB < 0 ? JROW >= J+2 : JROW <= J+2; JROW += -NNB) {
                     PPW = PPWO;
                     LEN  = 2 + J - JCOL;
                     for (I = JROW+NNB-1; I >= JROW; I--) {
                        CTEMP = A( I, J );
                        A[I][J] = CZERO;
                        S = B( I, J );
                        B[I][J] = CZERO;
                        for (JJ = PPW; JJ <= PPW+LEN-1; JJ++) {
                           TEMP = WORK( JJ + 2*NNB );
                           WORK[JJ + 2*NNB] = CTEMP*TEMP - CONJG( S )*WORK( JJ );
                           WORK[JJ] = S*TEMP + CTEMP*WORK( JJ );
                        }
                        LEN = LEN + 1;
                        PPW = PPW - 2*NNB - 1;
                     }
                     PPWO = PPWO + 4*NNB*NNB;
                  }
               }
            } else {

               claset('Lower', IHI - JCOL - 1, NNB, CZERO, CZERO, A( JCOL + 2, JCOL ), LDA );
               claset('Lower', IHI - JCOL - 1, NNB, CZERO, CZERO, B( JCOL + 2, JCOL ), LDB );
            }

            // Apply accumulated unitary matrices to A and B.

            if ( TOP > 0 ) {
               J = IHI - NBLST + 1;
               cgemm('No Transpose', 'No Transpose', TOP, NBLST, NBLST, CONE, A( 1, J ), LDA, WORK, NBLST, CZERO, WORK( PW ), TOP );
               clacpy('All', TOP, NBLST, WORK( PW ), TOP, A( 1, J ), LDA );
               PPWO = NBLST*NBLST + 1;
               J0 = J - NNB;
               for (J = J0; -NNB < 0 ? J >= JCOL+1 : J <= JCOL+1; J += -NNB) {
                  if ( BLK22 ) {

                     // Exploit the structure of U.

                     cunm22('Right', 'No Transpose', TOP, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, A( 1, J ), LDA, WORK( PW ), LWORK-PW+1, IERR );
                  } else {

                     // Ignore the structure of U.

                     cgemm('No Transpose', 'No Transpose', TOP, 2*NNB, 2*NNB, CONE, A( 1, J ), LDA, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), TOP );
                     clacpy('All', TOP, 2*NNB, WORK( PW ), TOP, A( 1, J ), LDA );
                  }
                  PPWO = PPWO + 4*NNB*NNB;
               }

               J = IHI - NBLST + 1;
               cgemm('No Transpose', 'No Transpose', TOP, NBLST, NBLST, CONE, B( 1, J ), LDB, WORK, NBLST, CZERO, WORK( PW ), TOP );
               clacpy('All', TOP, NBLST, WORK( PW ), TOP, B( 1, J ), LDB );
               PPWO = NBLST*NBLST + 1;
               J0 = J - NNB;
               for (J = J0; -NNB < 0 ? J >= JCOL+1 : J <= JCOL+1; J += -NNB) {
                  if ( BLK22 ) {

                     // Exploit the structure of U.

                     cunm22('Right', 'No Transpose', TOP, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, B( 1, J ), LDB, WORK( PW ), LWORK-PW+1, IERR );
                  } else {

                     // Ignore the structure of U.

                     cgemm('No Transpose', 'No Transpose', TOP, 2*NNB, 2*NNB, CONE, B( 1, J ), LDB, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), TOP );
                     clacpy('All', TOP, 2*NNB, WORK( PW ), TOP, B( 1, J ), LDB );
                  }
                  PPWO = PPWO + 4*NNB*NNB;
               }
            }

            // Apply accumulated unitary matrices to Z.

            if ( WANTZ ) {
               J = IHI - NBLST + 1;
               if ( INITQ ) {
                  TOPQ = max( 2, J - JCOL + 1 );
                  NH  = IHI - TOPQ + 1;
               } else {
                  TOPQ = 1;
                  NH = N;
               }
               cgemm('No Transpose', 'No Transpose', NH, NBLST, NBLST, CONE, Z( TOPQ, J ), LDZ, WORK, NBLST, CZERO, WORK( PW ), NH );
               clacpy('All', NH, NBLST, WORK( PW ), NH, Z( TOPQ, J ), LDZ );
               PPWO = NBLST*NBLST + 1;
               J0 = J - NNB;
               for (J = J0; -NNB < 0 ? J >= JCOL+1 : J <= JCOL+1; J += -NNB) {
                     if ( INITQ ) {
                     TOPQ = max( 2, J - JCOL + 1 );
                     NH  = IHI - TOPQ + 1;
                  }
                  if ( BLK22 ) {

                     // Exploit the structure of U.

                     cunm22('Right', 'No Transpose', NH, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, Z( TOPQ, J ), LDZ, WORK( PW ), LWORK-PW+1, IERR );
                  } else {

                     // Ignore the structure of U.

                     cgemm('No Transpose', 'No Transpose', NH, 2*NNB, 2*NNB, CONE, Z( TOPQ, J ), LDZ, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), NH );
                     clacpy('All', NH, 2*NNB, WORK( PW ), NH, Z( TOPQ, J ), LDZ );
                  }
                  PPWO = PPWO + 4*NNB*NNB;
               }
            }
         }
      }

      // Use unblocked code to reduce the rest of the matrix
      // Avoid re-initialization of modified Q and Z.

      COMPQ2 = COMPQ;
      COMPZ2 = COMPZ;
      if ( JCOL != ILO ) {
         if (WANTQ) COMPQ2 = 'V';
         IF ( WANTZ ) COMPZ2 = 'V';
      }

      if (JCOL < IHI) cgghrd( COMPQ2, COMPZ2, N, JCOL, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IERR );

      WORK[1] = SROUNDUP_LWORK( LWKOPT );

      }
