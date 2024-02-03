      SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, VECT;
      int                INFO, KD, LDAB, LDQ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               INITQ, UPPER, WANTQ;
      int                I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J, J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1, KDM1, KDN, L, LAST, LEND, NQ, NR, NRT;
      double             ABST;
      COMPLEX*16         T, TEMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLAR2V, ZLARGV, ZLARTG, ZLARTV, ZLASET, ZROT, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INITQ = LSAME( VECT, 'V' )
      WANTQ = INITQ .OR. LSAME( VECT, 'U' )
      UPPER = LSAME( UPLO, 'U' )
      KD1 = KD + 1
      KDM1 = KD - 1
      INCX = LDAB - 1
      IQEND = 1

      INFO = 0
      if ( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) {
         INFO = -1
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KD.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KD1 ) {
         INFO = -6
      } else if ( LDQ.LT.MAX( 1, N ) .AND. WANTQ ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZHBTRD', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize Q to the unit matrix, if needed

      IF( INITQ ) CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )

      // Wherever possible, plane rotations are generated and applied in
      // vector operations of length NR over the index set J1:J2:KD1.

      // The real cosines and complex sines of the plane rotations are
      // stored in the arrays D and WORK.

      INCA = KD1*LDAB
      KDN = MIN( N-1, KD )
      if ( UPPER ) {

         if ( KD.GT.1 ) {

            // Reduce to complex Hermitian tridiagonal form, working with
           t // he upper triangle

            NR = 0
            J1 = KDN + 2
            J2 = 1

            AB( KD1, 1 ) = DBLE( AB( KD1, 1 ) )
            DO 90 I = 1, N - 2

               // Reduce i-th row of matrix to tridiagonal form

               DO 80 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN

                  if ( NR.GT.0 ) {

                     // generate plane rotations to annihilate nonzero
                     // elements which have been created outside the band

                     CALL ZLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ), KD1, D( J1 ), KD1 )

                     // apply rotations from the right


                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                     if ( NR.GE.2*KD-1 ) {
                        DO 10 L = 1, KD - 1
                           CALL ZLARTV( NR, AB( L+1, J1-1 ), INCA, AB( L, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 )
   10                   CONTINUE

                     } else {
                        JEND = J1 + ( NR-1 )*KD1
                        DO 20 JINC = J1, JEND, KD1
                           CALL ZROT( KDM1, AB( 2, JINC-1 ), 1, AB( 1, JINC ), 1, D( JINC ), WORK( JINC ) )
   20                   CONTINUE
                     }
                  }


                  if ( K.GT.2 ) {
                     if ( K.LE.N-I+1 ) {

                        // generate plane rotation to annihilate a(i,i+k-1)
                        // within the band

                        CALL ZLARTG( AB( KD-K+3, I+K-2 ), AB( KD-K+2, I+K-1 ), D( I+K-1 ), WORK( I+K-1 ), TEMP )
                        AB( KD-K+3, I+K-2 ) = TEMP

                        // apply rotation from the right

                        CALL ZROT( K-3, AB( KD-K+4, I+K-2 ), 1, AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ), WORK( I+K-1 ) )
                     }
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  }

                  // apply plane rotations from both sides to diagonal
                  // blocks

                  IF( NR.GT.0 ) CALL ZLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ), AB( KD, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 )

                  // apply plane rotations from the left

                  if ( NR.GT.0 ) {
                     CALL ZLACGV( NR, WORK( J1 ), KD1 )
                     if ( 2*KD-1.LT.NR ) {

                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                        DO 30 L = 1, KD - 1
                           if ( J2+L.GT.N ) {
                              NRT = NR - 1
                           } else {
                              NRT = NR
                           }
                           IF( NRT.GT.0 ) CALL ZLARTV( NRT, AB( KD-L, J1+L ), INCA, AB( KD-L+1, J1+L ), INCA, D( J1 ), WORK( J1 ), KD1 )
   30                   CONTINUE
                     } else {
                        J1END = J1 + KD1*( NR-2 )
                        if ( J1END.GE.J1 ) {
                           DO 40 JIN = J1, J1END, KD1
                              CALL ZROT( KD-1, AB( KD-1, JIN+1 ), INCX, AB( KD, JIN+1 ), INCX, D( JIN ), WORK( JIN ) )
   40                      CONTINUE
                        }
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 ) CALL ZROT( LEND, AB( KD-1, LAST+1 ), INCX, AB( KD, LAST+1 ), INCX, D( LAST ), WORK( LAST ) )
                     }
                  }

                  if ( WANTQ ) {

                     // accumulate product of plane rotations in Q

                     if ( INITQ ) {

                 t // ake advantage of the fact that Q was
                  // initially the Identity matrix

                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 ) IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 50 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL ZROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), 1, D( J ), DCONJG( WORK( J ) ) )
   50                   CONTINUE
                     } else {

                        DO 60 J = J1, J2, KD1
                           CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1, D( J ), DCONJG( WORK( J ) ) )
   60                   CONTINUE
                     }

                  }

                  if ( J2+KDN.GT.N ) {

                     // adjust J2 to keep within the bounds of the matrix

                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  }

                  DO 70 J = J1, J2, KD1

                     // create nonzero element a(j-1,j+kd) outside the band
                     // and store it in WORK

                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                     AB( 1, J+KD ) = D( J )*AB( 1, J+KD )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
         }

         if ( KD.GT.0 ) {

            // make off-diagonal elements real and copy them to E

            DO 100 I = 1, N - 1
               T = AB( KD, I+1 )
               ABST = ABS( T )
               AB( KD, I+1 ) = ABST
               E( I ) = ABST
               if ( ABST.NE.ZERO ) {
                  T = T / ABST
               } else {
                  T = CONE
               }
               IF( I.LT.N-1 ) AB( KD, I+2 ) = AB( KD, I+2 )*T
               if ( WANTQ ) {
                  CALL ZSCAL( N, DCONJG( T ), Q( 1, I+1 ), 1 )
               }
  100       CONTINUE
         } else {

            // set E to zero if original matrix was diagonal

            DO 110 I = 1, N - 1
               E( I ) = ZERO
  110       CONTINUE
         }

         // copy diagonal elements to D

         DO 120 I = 1, N
            D( I ) = DBLE( AB( KD1, I ) )
  120    CONTINUE

      } else {

         if ( KD.GT.1 ) {

            // Reduce to complex Hermitian tridiagonal form, working with
           t // he lower triangle

            NR = 0
            J1 = KDN + 2
            J2 = 1

            AB( 1, 1 ) = DBLE( AB( 1, 1 ) )
            DO 210 I = 1, N - 2

               // Reduce i-th column of matrix to tridiagonal form

               DO 200 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN

                  if ( NR.GT.0 ) {

                     // generate plane rotations to annihilate nonzero
                     // elements which have been created outside the band

                     CALL ZLARGV( NR, AB( KD1, J1-KD1 ), INCA, WORK( J1 ), KD1, D( J1 ), KD1 )

                     // apply plane rotations from one side


                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                     if ( NR.GT.2*KD-1 ) {
                        DO 130 L = 1, KD - 1
                           CALL ZLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA, AB( KD1-L+1, J1-KD1+L ), INCA, D( J1 ), WORK( J1 ), KD1 )
  130                   CONTINUE
                     } else {
                        JEND = J1 + KD1*( NR-1 )
                        DO 140 JINC = J1, JEND, KD1
                           CALL ZROT( KDM1, AB( KD, JINC-KD ), INCX, AB( KD1, JINC-KD ), INCX, D( JINC ), WORK( JINC ) )
  140                   CONTINUE
                     }

                  }

                  if ( K.GT.2 ) {
                     if ( K.LE.N-I+1 ) {

                        // generate plane rotation to annihilate a(i+k-1,i)
                        // within the band

                        CALL ZLARTG( AB( K-1, I ), AB( K, I ), D( I+K-1 ), WORK( I+K-1 ), TEMP )
                        AB( K-1, I ) = TEMP

                        // apply rotation from the left

                        CALL ZROT( K-3, AB( K-2, I+1 ), LDAB-1, AB( K-1, I+1 ), LDAB-1, D( I+K-1 ), WORK( I+K-1 ) )
                     }
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  }

                  // apply plane rotations from both sides to diagonal
                  // blocks

                  IF( NR.GT.0 ) CALL ZLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ), AB( 2, J1-1 ), INCA, D( J1 ), WORK( J1 ), KD1 )

                  // apply plane rotations from the right


                     // Dependent on the the number of diagonals either
                     // ZLARTV or ZROT is used

                  if ( NR.GT.0 ) {
                     CALL ZLACGV( NR, WORK( J1 ), KD1 )
                     if ( NR.GT.2*KD-1 ) {
                        DO 150 L = 1, KD - 1
                           if ( J2+L.GT.N ) {
                              NRT = NR - 1
                           } else {
                              NRT = NR
                           }
                           IF( NRT.GT.0 ) CALL ZLARTV( NRT, AB( L+2, J1-1 ), INCA, AB( L+1, J1 ), INCA, D( J1 ), WORK( J1 ), KD1 )
  150                   CONTINUE
                     } else {
                        J1END = J1 + KD1*( NR-2 )
                        if ( J1END.GE.J1 ) {
                           DO 160 J1INC = J1, J1END, KD1
                              CALL ZROT( KDM1, AB( 3, J1INC-1 ), 1, AB( 2, J1INC ), 1, D( J1INC ), WORK( J1INC ) )
  160                      CONTINUE
                        }
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 ) CALL ZROT( LEND, AB( 3, LAST-1 ), 1, AB( 2, LAST ), 1, D( LAST ), WORK( LAST ) )
                     }
                  }



                  if ( WANTQ ) {

                     // accumulate product of plane rotations in Q

                     if ( INITQ ) {

                 t // ake advantage of the fact that Q was
                  // initially the Identity matrix

                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 ) IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 170 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL ZROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ), 1, D( J ), WORK( J ) )
  170                   CONTINUE
                     } else {

                        DO 180 J = J1, J2, KD1
                           CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1, D( J ), WORK( J ) )
  180                   CONTINUE
                     }
                  }

                  if ( J2+KDN.GT.N ) {

                     // adjust J2 to keep within the bounds of the matrix

                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  }

                  DO 190 J = J1, J2, KD1

                     // create nonzero element a(j+kd,j-1) outside the
                     // band and store it in WORK

                     WORK( J+KD ) = WORK( J )*AB( KD1, J )
                     AB( KD1, J ) = D( J )*AB( KD1, J )
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
         }

         if ( KD.GT.0 ) {

            // make off-diagonal elements real and copy them to E

            DO 220 I = 1, N - 1
               T = AB( 2, I )
               ABST = ABS( T )
               AB( 2, I ) = ABST
               E( I ) = ABST
               if ( ABST.NE.ZERO ) {
                  T = T / ABST
               } else {
                  T = CONE
               }
               IF( I.LT.N-1 ) AB( 2, I+1 ) = AB( 2, I+1 )*T
               if ( WANTQ ) {
                  CALL ZSCAL( N, T, Q( 1, I+1 ), 1 )
               }
  220       CONTINUE
         } else {

            // set E to zero if original matrix was diagonal

            DO 230 I = 1, N - 1
               E( I ) = ZERO
  230       CONTINUE
         }

         // copy diagonal elements to D

         DO 240 I = 1, N
            D( I ) = DBLE( AB( 1, I ) )
  240    CONTINUE
      }

      RETURN

      // End of ZHBTRD

      }
