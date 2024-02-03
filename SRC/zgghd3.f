      SUBROUTINE ZGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             COMPQ, COMPZ;
      int                IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE, CZERO
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      bool               BLK22, INITQ, INITZ, LQUERY, WANTQ, WANTZ;
      String             COMPQ2, COMPZ2;
      int                COLA, I, IERR, J, J0, JCOL, JJ, JROW, K, KACC22, LEN, LWKOPT, N2NB, NB, NBLST, NBMIN, NH, NNB, NX, PPW, PPWO, PW, TOP, TOPQ;
      double             C;
      COMPLEX*16         C1, C2, CTEMP, S, S1, S2, TEMP, TEMP1, TEMP2, TEMP3
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGGHRD, ZLARTG, ZLASET, ZUNM22, ZROT, ZGEMM, ZGEMV, ZTRMV, ZLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      INFO = 0
      NB = ILAENV( 1, 'ZGGHD3', ' ', N, ILO, IHI, -1 )
      NH = IHI - ILO + 1
      IF( NH.LE.1 ) THEN
         LWKOPT = 1
      ELSE
         LWKOPT = 6*N*NB
      END IF
      WORK( 1 ) = DCMPLX( LWKOPT )
      INITQ = LSAME( COMPQ, 'I' )
      WANTQ = INITQ .OR. LSAME( COMPQ, 'V' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      LQUERY = ( LWORK.EQ.-1 )

      IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( WANTQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( WANTZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGHD3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Initialize Q and Z if desired.

      IF( INITQ ) CALL ZLASET( 'All', N, N, CZERO, CONE, Q, LDQ )       IF( INITZ ) CALL ZLASET( 'All', N, N, CZERO, CONE, Z, LDZ )

      // Zero out lower triangle of B.

      IF( N.GT.1 ) CALL ZLASET( 'Lower', N-1, N-1, CZERO, CZERO, B(2, 1), LDB )

      // Quick return if possible

      IF( NH.LE.1 ) THEN
         WORK( 1 ) = CONE
         RETURN
      END IF

      // Determine the blocksize.

      NBMIN = ILAENV( 2, 'ZGGHD3', ' ', N, ILO, IHI, -1 )
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN

         // Determine when to use unblocked instead of blocked code.

         NX = MAX( NB, ILAENV( 3, 'ZGGHD3', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN

            // Determine if workspace is large enough for blocked code.

            IF( LWORK.LT.LWKOPT ) THEN

               // Not enough workspace to use optimal NB:  determine the
               // minimum value of NB, and reduce NB or force use of
               // unblocked code.

               NBMIN = MAX( 2, ILAENV( 2, 'ZGGHD3', ' ', N, ILO, IHI, -1 ) )
               IF( LWORK.GE.6*N*NBMIN ) THEN
                  NB = LWORK / ( 6*N )
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF

      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN

         // Use unblocked code below

         JCOL = ILO

      ELSE

         // Use blocked code

         KACC22 = ILAENV( 16, 'ZGGHD3', ' ', N, ILO, IHI, -1 )
         BLK22 = KACC22.EQ.2
         DO JCOL = ILO, IHI-2, NB
            NNB = MIN( NB, IHI-JCOL-1 )

            // Initialize small unitary factors that will hold the
            // accumulated Givens rotations in workspace.
            // N2NB   denotes the number of 2*NNB-by-2*NNB factors
            // NBLST  denotes the (possibly smaller) order of the last
                   // factor.

            N2NB = ( IHI-JCOL-1 ) / NNB - 1
            NBLST = IHI - JCOL - N2NB*NNB
            CALL ZLASET( 'All', NBLST, NBLST, CZERO, CONE, WORK, NBLST )
            PW = NBLST * NBLST + 1
            DO I = 1, N2NB
               CALL ZLASET( 'All', 2*NNB, 2*NNB, CZERO, CONE, WORK( PW ), 2*NNB )
               PW = PW + 4*NNB*NNB
            END DO

            // Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form.

            DO J = JCOL, JCOL+NNB-1

               // Reduce Jth column of A. Store cosines and sines in Jth
               // column of A and B, respectively.

               DO I = IHI, J+2, -1
                  TEMP = A( I-1, J )
                  CALL ZLARTG( TEMP, A( I, J ), C, S, A( I-1, J ) )
                  A( I, J ) = DCMPLX( C )
                  B( I, J ) = S
               END DO

               // Accumulate Givens rotations into workspace array.

               PPW  = ( NBLST + 1 )*( NBLST - 2 ) - J + JCOL + 1
               LEN  = 2 + J - JCOL
               JROW = J + N2NB*NNB + 2
               DO I = IHI, JROW, -1
                  CTEMP = A( I, J )
                  S = B( I, J )
                  DO JJ = PPW, PPW+LEN-1
                     TEMP = WORK( JJ + NBLST )
                     WORK( JJ + NBLST ) = CTEMP*TEMP - S*WORK( JJ )
                     WORK( JJ ) = DCONJG( S )*TEMP + CTEMP*WORK( JJ )
                  END DO
                  LEN = LEN + 1
                  PPW = PPW - NBLST - 1
               END DO

               PPWO = NBLST*NBLST + ( NNB+J-JCOL-1 )*2*NNB + NNB
               J0 = JROW - NNB
               DO JROW = J0, J+2, -NNB
                  PPW = PPWO
                  LEN  = 2 + J - JCOL
                  DO I = JROW+NNB-1, JROW, -1
                     CTEMP = A( I, J )
                     S = B( I, J )
                     DO JJ = PPW, PPW+LEN-1
                        TEMP = WORK( JJ + 2*NNB )
                        WORK( JJ + 2*NNB ) = CTEMP*TEMP - S*WORK( JJ )
                        WORK( JJ ) = DCONJG( S )*TEMP + CTEMP*WORK( JJ )
                     END DO
                     LEN = LEN + 1
                     PPW = PPW - 2*NNB - 1
                  END DO
                  PPWO = PPWO + 4*NNB*NNB
               END DO

               // TOP denotes the number of top rows in A and B that will
               // not be updated during the next steps.

               IF( JCOL.LE.2 ) THEN
                  TOP = 0
               ELSE
                  TOP = JCOL
               END IF

               // Propagate transformations through B and replace stored
               // left sines/cosines by right sines/cosines.

               DO JJ = N, J+1, -1

                  // Update JJth column of B.

                  DO I = MIN( JJ+1, IHI ), J+2, -1
                     CTEMP = A( I, J )
                     S = B( I, J )
                     TEMP = B( I, JJ )
                     B( I, JJ ) = CTEMP*TEMP - DCONJG( S )*B( I-1, JJ )
                     B( I-1, JJ ) = S*TEMP + CTEMP*B( I-1, JJ )
                  END DO

                  // Annihilate B( JJ+1, JJ ).

                  IF( JJ.LT.IHI ) THEN
                     TEMP = B( JJ+1, JJ+1 )
                     CALL ZLARTG( TEMP, B( JJ+1, JJ ), C, S, B( JJ+1, JJ+1 ) )
                     B( JJ+1, JJ ) = CZERO
                     CALL ZROT( JJ-TOP, B( TOP+1, JJ+1 ), 1, B( TOP+1, JJ ), 1, C, S )
                     A( JJ+1, J ) = DCMPLX( C )
                     B( JJ+1, J ) = -DCONJG( S )
                  END IF
               END DO

               // Update A by transformations from right.

               JJ = MOD( IHI-J-1, 3 )
               DO I = IHI-J-3, JJ+1, -3
                  CTEMP = A( J+1+I, J )
                  S = -B( J+1+I, J )
                  C1 = A( J+2+I, J )
                  S1 = -B( J+2+I, J )
                  C2 = A( J+3+I, J )
                  S2 = -B( J+3+I, J )

                  DO K = TOP+1, IHI
                     TEMP = A( K, J+I  )
                     TEMP1 = A( K, J+I+1 )
                     TEMP2 = A( K, J+I+2 )
                     TEMP3 = A( K, J+I+3 )
                     A( K, J+I+3 ) = C2*TEMP3 + DCONJG( S2 )*TEMP2
                     TEMP2 = -S2*TEMP3 + C2*TEMP2
                     A( K, J+I+2 ) = C1*TEMP2 + DCONJG( S1 )*TEMP1
                     TEMP1 = -S1*TEMP2 + C1*TEMP1
                     A( K, J+I+1 ) = CTEMP*TEMP1 + DCONJG( S )*TEMP
                     A( K, J+I ) = -S*TEMP1 + CTEMP*TEMP
                  END DO
               END DO

               IF( JJ.GT.0 ) THEN
                  DO I = JJ, 1, -1
                     C = DBLE( A( J+1+I, J ) )
                     CALL ZROT( IHI-TOP, A( TOP+1, J+I+1 ), 1, A( TOP+1, J+I ), 1, C, -DCONJG( B( J+1+I, J ) ) )
                  END DO
               END IF

               // Update (J+1)th column of A by transformations from left.

               IF ( J .LT. JCOL + NNB - 1 ) THEN
                  LEN  = 1 + J - JCOL

                  // Multiply with the trailing accumulated unitary
                  // matrix, which takes the form

                         // [  U11  U12  ]
                     // U = [            ],
                         // [  U21  U22  ]

                  // where U21 is a LEN-by-LEN matrix and U12 is lower
                 t // riangular.

                  JROW = IHI - NBLST + 1
                  CALL ZGEMV( 'Conjugate', NBLST, LEN, CONE, WORK, NBLST, A( JROW, J+1 ), 1, CZERO, WORK( PW ), 1 )
                  PPW = PW + LEN
                  DO I = JROW, JROW+NBLST-LEN-1
                     WORK( PPW ) = A( I, J+1 )
                     PPW = PPW + 1
                  END DO
                  CALL ZTRMV( 'Lower', 'Conjugate', 'Non-unit', NBLST-LEN, WORK( LEN*NBLST + 1 ), NBLST, WORK( PW+LEN ), 1 )                   CALL ZGEMV( 'Conjugate', LEN, NBLST-LEN, CONE, WORK( (LEN+1)*NBLST - LEN + 1 ), NBLST, A( JROW+NBLST-LEN, J+1 ), 1, CONE, WORK( PW+LEN ), 1 )
                  PPW = PW
                  DO I = JROW, JROW+NBLST-1
                     A( I, J+1 ) = WORK( PPW )
                     PPW = PPW + 1
                  END DO

                  // Multiply with the other accumulated unitary
                  // matrices, which take the form

                         // [  U11  U12   0  ]
                         // [                ]
                     // U = [  U21  U22   0  ],
                         // [                ]
                         // [   0    0    I  ]

                  // where I denotes the (NNB-LEN)-by-(NNB-LEN) identity
                  // matrix, U21 is a LEN-by-LEN upper triangular matrix
                  // and U12 is an NNB-by-NNB lower triangular matrix.

                  PPWO = 1 + NBLST*NBLST
                  J0 = JROW - NNB
                  DO JROW = J0, JCOL+1, -NNB
                     PPW = PW + LEN
                     DO I = JROW, JROW+NNB-1
                        WORK( PPW ) = A( I, J+1 )
                        PPW = PPW + 1
                     END DO
                     PPW = PW
                     DO I = JROW+NNB, JROW+NNB+LEN-1
                        WORK( PPW ) = A( I, J+1 )
                        PPW = PPW + 1
                     END DO
                     CALL ZTRMV( 'Upper', 'Conjugate', 'Non-unit', LEN, WORK( PPWO + NNB ), 2*NNB, WORK( PW ), 1 )                      CALL ZTRMV( 'Lower', 'Conjugate', 'Non-unit', NNB, WORK( PPWO + 2*LEN*NNB ), 2*NNB, WORK( PW + LEN ), 1 )                      CALL ZGEMV( 'Conjugate', NNB, LEN, CONE, WORK( PPWO ), 2*NNB, A( JROW, J+1 ), 1, CONE, WORK( PW ), 1 )                      CALL ZGEMV( 'Conjugate', LEN, NNB, CONE, WORK( PPWO + 2*LEN*NNB + NNB ), 2*NNB, A( JROW+NNB, J+1 ), 1, CONE, WORK( PW+LEN ), 1 )
                     PPW = PW
                     DO I = JROW, JROW+LEN+NNB-1
                        A( I, J+1 ) = WORK( PPW )
                        PPW = PPW + 1
                     END DO
                     PPWO = PPWO + 4*NNB*NNB
                  END DO
               END IF
            END DO

            // Apply accumulated unitary matrices to A.

            COLA = N - JCOL - NNB + 1
            J = IHI - NBLST + 1
            CALL ZGEMM( 'Conjugate', 'No Transpose', NBLST, COLA, NBLST, CONE, WORK, NBLST, A( J, JCOL+NNB ), LDA, CZERO, WORK( PW ), NBLST )
            CALL ZLACPY( 'All', NBLST, COLA, WORK( PW ), NBLST, A( J, JCOL+NNB ), LDA )
            PPWO = NBLST*NBLST + 1
            J0 = J - NNB
            DO J = J0, JCOL+1, -NNB
               IF ( BLK22 ) THEN

                  // Exploit the structure of

                         // [  U11  U12  ]
                     // U = [            ]
                         // [  U21  U22  ],

                  // where all blocks are NNB-by-NNB, U21 is upper
                 t // riangular and U12 is lower triangular.

                  CALL ZUNM22( 'Left', 'Conjugate', 2*NNB, COLA, NNB, NNB, WORK( PPWO ), 2*NNB, A( J, JCOL+NNB ), LDA, WORK( PW ), LWORK-PW+1, IERR )
               ELSE

                  // Ignore the structure of U.

                  CALL ZGEMM( 'Conjugate', 'No Transpose', 2*NNB, COLA, 2*NNB, CONE, WORK( PPWO ), 2*NNB, A( J, JCOL+NNB ), LDA, CZERO, WORK( PW ), 2*NNB )
                  CALL ZLACPY( 'All', 2*NNB, COLA, WORK( PW ), 2*NNB, A( J, JCOL+NNB ), LDA )
               END IF
               PPWO = PPWO + 4*NNB*NNB
            END DO

            // Apply accumulated unitary matrices to Q.

            IF( WANTQ ) THEN
               J = IHI - NBLST + 1
               IF ( INITQ ) THEN
                  TOPQ = MAX( 2, J - JCOL + 1 )
                  NH  = IHI - TOPQ + 1
               ELSE
                  TOPQ = 1
                  NH = N
               END IF
               CALL ZGEMM( 'No Transpose', 'No Transpose', NH, NBLST, NBLST, CONE, Q( TOPQ, J ), LDQ, WORK, NBLST, CZERO, WORK( PW ), NH )
               CALL ZLACPY( 'All', NH, NBLST, WORK( PW ), NH, Q( TOPQ, J ), LDQ )
               PPWO = NBLST*NBLST + 1
               J0 = J - NNB
               DO J = J0, JCOL+1, -NNB
                  IF ( INITQ ) THEN
                     TOPQ = MAX( 2, J - JCOL + 1 )
                     NH  = IHI - TOPQ + 1
                  END IF
                  IF ( BLK22 ) THEN

                     // Exploit the structure of U.

                     CALL ZUNM22( 'Right', 'No Transpose', NH, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, Q( TOPQ, J ), LDQ, WORK( PW ), LWORK-PW+1, IERR )
                  ELSE

                     // Ignore the structure of U.

                     CALL ZGEMM( 'No Transpose', 'No Transpose', NH, 2*NNB, 2*NNB, CONE, Q( TOPQ, J ), LDQ, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), NH )
                     CALL ZLACPY( 'All', NH, 2*NNB, WORK( PW ), NH, Q( TOPQ, J ), LDQ )
                  END IF
                  PPWO = PPWO + 4*NNB*NNB
               END DO
            END IF

            // Accumulate right Givens rotations if required.

            IF ( WANTZ .OR. TOP.GT.0 ) THEN

               // Initialize small unitary factors that will hold the
               // accumulated Givens rotations in workspace.

               CALL ZLASET( 'All', NBLST, NBLST, CZERO, CONE, WORK, NBLST )
               PW = NBLST * NBLST + 1
               DO I = 1, N2NB
                  CALL ZLASET( 'All', 2*NNB, 2*NNB, CZERO, CONE, WORK( PW ), 2*NNB )
                  PW = PW + 4*NNB*NNB
               END DO

               // Accumulate Givens rotations into workspace array.

               DO J = JCOL, JCOL+NNB-1
                  PPW  = ( NBLST + 1 )*( NBLST - 2 ) - J + JCOL + 1
                  LEN  = 2 + J - JCOL
                  JROW = J + N2NB*NNB + 2
                  DO I = IHI, JROW, -1
                     CTEMP = A( I, J )
                     A( I, J ) = CZERO
                     S = B( I, J )
                     B( I, J ) = CZERO
                     DO JJ = PPW, PPW+LEN-1
                        TEMP = WORK( JJ + NBLST )
                        WORK( JJ + NBLST ) = CTEMP*TEMP - DCONJG( S )*WORK( JJ )
                        WORK( JJ ) = S*TEMP + CTEMP*WORK( JJ )
                     END DO
                     LEN = LEN + 1
                     PPW = PPW - NBLST - 1
                  END DO

                  PPWO = NBLST*NBLST + ( NNB+J-JCOL-1 )*2*NNB + NNB
                  J0 = JROW - NNB
                  DO JROW = J0, J+2, -NNB
                     PPW = PPWO
                     LEN  = 2 + J - JCOL
                     DO I = JROW+NNB-1, JROW, -1
                        CTEMP = A( I, J )
                        A( I, J ) = CZERO
                        S = B( I, J )
                        B( I, J ) = CZERO
                        DO JJ = PPW, PPW+LEN-1
                           TEMP = WORK( JJ + 2*NNB )
                           WORK( JJ + 2*NNB ) = CTEMP*TEMP - DCONJG( S )*WORK( JJ )
                           WORK( JJ ) = S*TEMP + CTEMP*WORK( JJ )
                        END DO
                        LEN = LEN + 1
                        PPW = PPW - 2*NNB - 1
                     END DO
                     PPWO = PPWO + 4*NNB*NNB
                  END DO
               END DO
            ELSE

               CALL ZLASET( 'Lower', IHI - JCOL - 1, NNB, CZERO, CZERO, A( JCOL + 2, JCOL ), LDA )                CALL ZLASET( 'Lower', IHI - JCOL - 1, NNB, CZERO, CZERO, B( JCOL + 2, JCOL ), LDB )
            END IF

            // Apply accumulated unitary matrices to A and B.

            IF ( TOP.GT.0 ) THEN
               J = IHI - NBLST + 1
               CALL ZGEMM( 'No Transpose', 'No Transpose', TOP, NBLST, NBLST, CONE, A( 1, J ), LDA, WORK, NBLST, CZERO, WORK( PW ), TOP )
               CALL ZLACPY( 'All', TOP, NBLST, WORK( PW ), TOP, A( 1, J ), LDA )
               PPWO = NBLST*NBLST + 1
               J0 = J - NNB
               DO J = J0, JCOL+1, -NNB
                  IF ( BLK22 ) THEN

                     // Exploit the structure of U.

                     CALL ZUNM22( 'Right', 'No Transpose', TOP, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, A( 1, J ), LDA, WORK( PW ), LWORK-PW+1, IERR )
                  ELSE

                     // Ignore the structure of U.

                     CALL ZGEMM( 'No Transpose', 'No Transpose', TOP, 2*NNB, 2*NNB, CONE, A( 1, J ), LDA, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), TOP )
                     CALL ZLACPY( 'All', TOP, 2*NNB, WORK( PW ), TOP, A( 1, J ), LDA )
                  END IF
                  PPWO = PPWO + 4*NNB*NNB
               END DO

               J = IHI - NBLST + 1
               CALL ZGEMM( 'No Transpose', 'No Transpose', TOP, NBLST, NBLST, CONE, B( 1, J ), LDB, WORK, NBLST, CZERO, WORK( PW ), TOP )
               CALL ZLACPY( 'All', TOP, NBLST, WORK( PW ), TOP, B( 1, J ), LDB )
               PPWO = NBLST*NBLST + 1
               J0 = J - NNB
               DO J = J0, JCOL+1, -NNB
                  IF ( BLK22 ) THEN

                     // Exploit the structure of U.

                     CALL ZUNM22( 'Right', 'No Transpose', TOP, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, B( 1, J ), LDB, WORK( PW ), LWORK-PW+1, IERR )
                  ELSE

                     // Ignore the structure of U.

                     CALL ZGEMM( 'No Transpose', 'No Transpose', TOP, 2*NNB, 2*NNB, CONE, B( 1, J ), LDB, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), TOP )
                     CALL ZLACPY( 'All', TOP, 2*NNB, WORK( PW ), TOP, B( 1, J ), LDB )
                  END IF
                  PPWO = PPWO + 4*NNB*NNB
               END DO
            END IF

            // Apply accumulated unitary matrices to Z.

            IF( WANTZ ) THEN
               J = IHI - NBLST + 1
               IF ( INITQ ) THEN
                  TOPQ = MAX( 2, J - JCOL + 1 )
                  NH  = IHI - TOPQ + 1
               ELSE
                  TOPQ = 1
                  NH = N
               END IF
               CALL ZGEMM( 'No Transpose', 'No Transpose', NH, NBLST, NBLST, CONE, Z( TOPQ, J ), LDZ, WORK, NBLST, CZERO, WORK( PW ), NH )
               CALL ZLACPY( 'All', NH, NBLST, WORK( PW ), NH, Z( TOPQ, J ), LDZ )
               PPWO = NBLST*NBLST + 1
               J0 = J - NNB
               DO J = J0, JCOL+1, -NNB
                     IF ( INITQ ) THEN
                     TOPQ = MAX( 2, J - JCOL + 1 )
                     NH  = IHI - TOPQ + 1
                  END IF
                  IF ( BLK22 ) THEN

                     // Exploit the structure of U.

                     CALL ZUNM22( 'Right', 'No Transpose', NH, 2*NNB, NNB, NNB, WORK( PPWO ), 2*NNB, Z( TOPQ, J ), LDZ, WORK( PW ), LWORK-PW+1, IERR )
                  ELSE

                     // Ignore the structure of U.

                     CALL ZGEMM( 'No Transpose', 'No Transpose', NH, 2*NNB, 2*NNB, CONE, Z( TOPQ, J ), LDZ, WORK( PPWO ), 2*NNB, CZERO, WORK( PW ), NH )
                     CALL ZLACPY( 'All', NH, 2*NNB, WORK( PW ), NH, Z( TOPQ, J ), LDZ )
                  END IF
                  PPWO = PPWO + 4*NNB*NNB
               END DO
            END IF
         END DO
      END IF

      // Use unblocked code to reduce the rest of the matrix
      // Avoid re-initialization of modified Q and Z.

      COMPQ2 = COMPQ
      COMPZ2 = COMPZ
      IF ( JCOL.NE.ILO ) THEN
         IF ( WANTQ ) COMPQ2 = 'V'          IF ( WANTZ ) COMPZ2 = 'V'
      END IF

      IF ( JCOL.LT.IHI ) CALL ZGGHRD( COMPQ2, COMPZ2, N, JCOL, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IERR )

      WORK( 1 ) = DCMPLX( LWKOPT )

      RETURN

      // End of ZGGHD3

      END
