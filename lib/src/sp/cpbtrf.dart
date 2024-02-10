      void cpbtrf(UPLO, N, KD, AB, LDAB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, N;
      Complex            AB( LDAB, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                NBMAX, LDWORK;
      const              NBMAX = 32, LDWORK = NBMAX+1 ;
      int                I, I2, I3, IB, II, J, JJ, NB;
      Complex            WORK( LDWORK, NBMAX );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CPBTF2, CPOTF2, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      // Test the input parameters.

      INFO = 0;
      if ( ( !lsame( UPLO, 'U' ) ) && ( !lsame( UPLO, 'L' ) ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( LDAB < KD+1 ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('CPBTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment

      NB = ilaenv( 1, 'CPBTRF', UPLO, N, KD, -1, -1 );

      // The block size must not exceed the semi-bandwidth KD, and must not
      // exceed the limit set by the size of the local array WORK.

      NB = min( NB, NBMAX );

      if ( NB <= 1 || NB > KD ) {

         // Use unblocked code

         cpbtf2(UPLO, N, KD, AB, LDAB, INFO );
      } else {

         // Use blocked code

         if ( lsame( UPLO, 'U' ) ) {

            // Compute the Cholesky factorization of a Hermitian band
            // matrix, given the upper triangle of the matrix in band
            // storage.

            // Zero the upper triangle of the work array.

            for (J = 1; J <= NB; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  WORK[I][J] = ZERO;
               } // 10
            } // 20

            // Process the band matrix one diagonal block at a time.

            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) { // 70
               IB = min( NB, N-I+1 );

               // Factorize the diagonal block

               cpotf2(UPLO, IB, AB( KD+1, I ), LDAB-1, II );
               if ( II != 0 ) {
                  INFO = I + II - 1;
                  GO TO 150;
               }
               if ( I+IB <= N ) {

                  // Update the relevant part of the trailing submatrix.
                  // If A11 denotes the diagonal block which has just been
                  // factorized, then we need to update the remaining
                  // blocks in the diagram:

                     // A11   A12   A13
                     //       A22   A23
                     //             A33

                  // The numbers of rows and columns in the partitioning
                  // are IB, I2, I3 respectively. The blocks A12, A22 and
                  // A23 are empty if IB = KD. The upper triangle of A13
                  // lies outside the band.

                  I2 = min( KD-IB, N-I-IB+1 );
                  I3 = min( IB, N-I-KD+1 );

                  if ( I2 > 0 ) {

                     // Update A12

                     ctrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', IB, I2, CONE, AB( KD+1, I ), LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 );

                     // Update A22

                     cherk('Upper', 'Conjugate transpose', I2, IB, -ONE, AB( KD+1-IB, I+IB ), LDAB-1, ONE, AB( KD+1, I+IB ), LDAB-1 );
                  }

                  if ( I3 > 0 ) {

                     // Copy the lower triangle of A13 into the work array.

                     for (JJ = 1; JJ <= I3; JJ++) { // 40
                        for (II = JJ; II <= IB; II++) { // 30
                           WORK[II][JJ] = AB( II-JJ+1, JJ+I+KD-1 );
                        } // 30
                     } // 40

                     // Update A13 (in the work array).

                     ctrsm('Left', 'Upper', 'Conjugate transpose', 'Non-unit', IB, I3, CONE, AB( KD+1, I ), LDAB-1, WORK, LDWORK );

                     // Update A23

                     if (I2 > 0) cgemm( 'Conjugate transpose', 'No transpose', I2, I3, IB, -CONE, AB( KD+1-IB, I+IB ), LDAB-1, WORK, LDWORK, CONE, AB( 1+IB, I+KD ), LDAB-1 );

                     // Update A33

                     cherk('Upper', 'Conjugate transpose', I3, IB, -ONE, WORK, LDWORK, ONE, AB( KD+1, I+KD ), LDAB-1 );

                     // Copy the lower triangle of A13 back into place.

                     for (JJ = 1; JJ <= I3; JJ++) { // 60
                        for (II = JJ; II <= IB; II++) { // 50
                           AB[II-JJ+1][JJ+I+KD-1] = WORK( II, JJ );
                        } // 50
                     } // 60
                  }
               }
            } // 70
         } else {

            // Compute the Cholesky factorization of a Hermitian band
            // matrix, given the lower triangle of the matrix in band
            // storage.

            // Zero the lower triangle of the work array.

            for (J = 1; J <= NB; J++) { // 90
               for (I = J + 1; I <= NB; I++) { // 80
                  WORK[I][J] = ZERO;
               } // 80
            } // 90

            // Process the band matrix one diagonal block at a time.

            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) { // 140
               IB = min( NB, N-I+1 );

               // Factorize the diagonal block

               cpotf2(UPLO, IB, AB( 1, I ), LDAB-1, II );
               if ( II != 0 ) {
                  INFO = I + II - 1;
                  GO TO 150;
               }
               if ( I+IB <= N ) {

                  // Update the relevant part of the trailing submatrix.
                  // If A11 denotes the diagonal block which has just been
                  // factorized, then we need to update the remaining
                  // blocks in the diagram:

                     // A11
                     // A21   A22
                     // A31   A32   A33

                  // The numbers of rows and columns in the partitioning
                  // are IB, I2, I3 respectively. The blocks A21, A22 and
                  // A32 are empty if IB = KD. The lower triangle of A31
                  // lies outside the band.

                  I2 = min( KD-IB, N-I-IB+1 );
                  I3 = min( IB, N-I-KD+1 );

                  if ( I2 > 0 ) {

                     // Update A21

                     ctrsm('Right', 'Lower', 'Conjugate transpose', 'Non-unit', I2, IB, CONE, AB( 1, I ), LDAB-1, AB( 1+IB, I ), LDAB-1 );

                     // Update A22

                     cherk('Lower', 'No transpose', I2, IB, -ONE, AB( 1+IB, I ), LDAB-1, ONE, AB( 1, I+IB ), LDAB-1 );
                  }

                  if ( I3 > 0 ) {

                     // Copy the upper triangle of A31 into the work array.

                     for (JJ = 1; JJ <= IB; JJ++) { // 110
                        for (II = 1; II <= min( JJ, I3 ); II++) { // 100
                           WORK[II][JJ] = AB( KD+1-JJ+II, JJ+I-1 );
                        } // 100
                     } // 110

                     // Update A31 (in the work array).

                     ctrsm('Right', 'Lower', 'Conjugate transpose', 'Non-unit', I3, IB, CONE, AB( 1, I ), LDAB-1, WORK, LDWORK );

                     // Update A32

                     if (I2 > 0) cgemm( 'No transpose', 'Conjugate transpose', I3, I2, IB, -CONE, WORK, LDWORK, AB( 1+IB, I ), LDAB-1, CONE, AB( 1+KD-IB, I+IB ), LDAB-1 );

                     // Update A33

                     cherk('Lower', 'No transpose', I3, IB, -ONE, WORK, LDWORK, ONE, AB( 1, I+KD ), LDAB-1 );

                     // Copy the upper triangle of A31 back into place.

                     for (JJ = 1; JJ <= IB; JJ++) { // 130
                        for (II = 1; II <= min( JJ, I3 ); II++) { // 120
                           AB[KD+1-JJ+II][JJ+I-1] = WORK( II, JJ );
                        } // 120
                     } // 130
                  }
               }
            } // 140
         }
      }
      return;

      } // 150
      }
