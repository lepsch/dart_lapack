      SUBROUTINE SDRVST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1, WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK, IWORK, LIWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               A( LDA, * ), D1( * ), D2( * ), D3( * ), D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ), U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ), WA3( * ), WORK( * ), Z( LDU, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, TEN
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, TEN = 10.0E0 ;
      REAL               HALF
      const              HALF = 0.5E+0 ;
      int                MAXTYP;
      const              MAXTYP = 18 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IDIAG, IHBW, IINFO, IL, IMODE, INDX, IROW, ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL, JSIZE, JTYPE, KD, LGN, LIWEDC, LWEDC, M, M2, M3, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLARND, SSXT1
      // EXTERNAL SLAMCH, SLARND, SSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, SLACPY, SLAFTS, SLASET, SLATMR, SLATMS, SSBEV, SSBEVD, SSBEVX, SSPEV, SSPEVD, SSPEVX, SSTEV, SSTEVD, SSTEVR, SSTEVX, SSTT21, SSTT22, SSYEV, SSYEVD, SSYEVR, SSYEVX, SSYT21, SSYEVD_2STAGE, SSYEVR_2STAGE, SSYEVX_2STAGE, SSYEV_2STAGE, SSBEV_2STAGE, SSBEVD_2STAGE, SSBEVX_2STAGE, SSYTRD_2STAGE, SSYTRD_SY2SB, SSYTRD_SB2ST, SSYT22, XERBLA
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 3*9 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4 /
      // ..
      // .. Executable Statements ..

      // Keep ftrnchek happy

      VL = ZERO
      VU = ZERO

      // 1)      Check for errors

      NTESTT = 0
      INFO = 0

      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDU.LT.NMAX ) {
         INFO = -16
      } else if ( 2*MAX( 2, NMAX )**2.GT.LWORK ) {
         INFO = -21
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SDRVST2STG', -INFO )
         RETURN
      }

      // Quick return if nothing to do

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      // Loop over sizes, types

      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
         ISEED3( I ) = ISEED( I )
   20 CONTINUE

      NERRS = 0
      NMATS = 0


      DO 1740 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         if ( N.GT.0 ) {
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2
            // LIWEDC = 6 + 6*N + 5*N*LGN
            LIWEDC = 3 + 5*N
         } else {
            LWEDC = 9
            // LIWEDC = 12
            LIWEDC = 8
         }
         ANINV = ONE / REAL( MAX( 1, N ) )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 1730 JTYPE = 1, MTYPES

            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 1730
            NMATS = NMATS + 1
            NTEST = 0

            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE

            // 2)      Compute "A"

                    // Control parameters:

                // KMAGN  KMODE        KTYPE
            // =1  O(1)   clustered 1  zero
            // =2  large  clustered 2  identity
            // =3  small  exponential  (none)
            // =4         arithmetic   diagonal, (w/ eigenvalues)
            // =5         random log   symmetric, w/ eigenvalues
            // =6         random       (none)
            // =7                      random diagonal
            // =8                      random symmetric
            // =9                      band symmetric, w/ eigenvalues

            IF( MTYPES.GT.MAXTYP ) GO TO 110

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40       CONTINUE
            ANORM = ONE
            GO TO 70

   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70

   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70

   70       CONTINUE

            CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            COND = ULPINV

            // Special Matrices -- Identity & Jordan block

                    // Zero

            if ( ITYPE.EQ.1 ) {
               IINFO = 0

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO )

            } else if ( ITYPE.EQ.5 ) {

               // Symmetric, eigenvalues specified

               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO )

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               IDUMMA( 1 ) = 1
               CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.8 ) {

               // Symmetric, random eigenvalues

               IDUMMA( 1 ) = 1
               CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.9 ) {

               // Symmetric banded, eigenvalues specified

               IHBW = INT( ( N-1 )*SLARND( 1, ISEED3 ) )
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, IHBW, IHBW, 'Z', U, LDU, WORK( N+1 ), IINFO )

               // Store as dense matrix for most routines.

               CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               DO 100 IDIAG = -IHBW, IHBW
                  IROW = IHBW - IDIAG + 1
                  J1 = MAX( 1, IDIAG+1 )
                  J2 = MIN( N, N+IDIAG )
                  DO 90 J = J1, J2
                     I = J - IDIAG
                     A( I, J ) = U( IROW, J )
   90             CONTINUE
  100          CONTINUE
            } else {
               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

  110       CONTINUE

            ABSTOL = UNFL + UNFL
            if ( N.LE.1 ) {
               IL = 1
               IU = N
            } else {
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               if ( IL.GT.IU ) {
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               }
            }

            // 3)      If matrix is tridiagonal, call SSTEV and SSTEVX.

            if ( JTYPE.LE.7 ) {
               NTEST = 1
               DO 120 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  120          CONTINUE
               DO 130 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  130          CONTINUE
               SRNAMT = 'SSTEV'
               CALL SSTEV( 'V', N, D1, D2, Z, LDU, WORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEV(V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 1 ) = ULPINV
                     RESULT( 2 ) = ULPINV
                     RESULT( 3 ) = ULPINV
                     GO TO 180
                  }
               }

               // Do tests 1 and 2.

               DO 140 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  140          CONTINUE
               DO 150 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  150          CONTINUE
               CALL SSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK, RESULT( 1 ) )

               NTEST = 3
               DO 160 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  160          CONTINUE
               SRNAMT = 'SSTEV'
               CALL SSTEV( 'N', N, D3, D4, Z, LDU, WORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEV(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 3 ) = ULPINV
                     GO TO 180
                  }
               }

               // Do test 3.

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 170 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  170          CONTINUE
               RESULT( 3 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  180          CONTINUE

               NTEST = 4
               DO 190 I = 1, N
                  EVEIGS( I ) = D3( I )
                  D1( I ) = REAL( A( I, I ) )
  190          CONTINUE
               DO 200 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  200          CONTINUE
               SRNAMT = 'SSTEVX'
               CALL SSTEVX( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 4 ) = ULPINV
                     RESULT( 5 ) = ULPINV
                     RESULT( 6 ) = ULPINV
                     GO TO 250
                  }
               }
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }

               // Do tests 4 and 5.

               DO 210 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  210          CONTINUE
               DO 220 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  220          CONTINUE
               CALL SSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK, RESULT( 4 ) )

               NTEST = 6
               DO 230 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  230          CONTINUE
               SRNAMT = 'SSTEVX'
               CALL SSTEVX( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 6 ) = ULPINV
                     GO TO 250
                  }
               }

               // Do test 6.

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 240 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), ABS( EVEIGS( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) )
  240          CONTINUE
               RESULT( 6 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  250          CONTINUE

               NTEST = 7
               DO 260 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  260          CONTINUE
               DO 270 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  270          CONTINUE
               SRNAMT = 'SSTEVR'
               CALL SSTEVR( 'V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 7 ) = ULPINV
                     RESULT( 8 ) = ULPINV
                     GO TO 320
                  }
               }
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }

               // Do tests 7 and 8.

               DO 280 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  280          CONTINUE
               DO 290 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  290          CONTINUE
               CALL SSTT21( N, 0, D3, D4, WA1, D2, Z, LDU, WORK, RESULT( 7 ) )

               NTEST = 9
               DO 300 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  300          CONTINUE
               SRNAMT = 'SSTEVR'
               CALL SSTEVR( 'N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 9 ) = ULPINV
                     GO TO 320
                  }
               }

               // Do test 9.

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 310 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), ABS( EVEIGS( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) )
  310          CONTINUE
               RESULT( 9 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  320          CONTINUE


               NTEST = 10
               DO 330 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  330          CONTINUE
               DO 340 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  340          CONTINUE
               SRNAMT = 'SSTEVX'
               CALL SSTEVX( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,I)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 10 ) = ULPINV
                     RESULT( 11 ) = ULPINV
                     RESULT( 12 ) = ULPINV
                     GO TO 380
                  }
               }

               // Do tests 10 and 11.

               DO 350 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  350          CONTINUE
               DO 360 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  360          CONTINUE
               CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, MAX( 1, M2 ), RESULT( 10 ) )


               NTEST = 12
               DO 370 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  370          CONTINUE
               SRNAMT = 'SSTEVX'
               CALL SSTEVX( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,I)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 12 ) = ULPINV
                     GO TO 380
                  }
               }

               // Do test 12.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 12 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )

  380          CONTINUE

               NTEST = 12
               if ( N.GT.0 ) {
                  if ( IL.NE.1 ) {
                     VL = WA1( IL ) - MAX( HALF* ( WA1( IL )-WA1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else {
                     VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
                  if ( IU.NE.N ) {
                     VU = WA1( IU ) + MAX( HALF* ( WA1( IU+1 )-WA1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else {
                     VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
               } else {
                  VL = ZERO
                  VU = ONE
               }

               DO 390 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  390          CONTINUE
               DO 400 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  400          CONTINUE
               SRNAMT = 'SSTEVX'
               CALL SSTEVX( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 13 ) = ULPINV
                     RESULT( 14 ) = ULPINV
                     RESULT( 15 ) = ULPINV
                     GO TO 440
                  }
               }

               if ( M2.EQ.0 .AND. N.GT.0 ) {
                  RESULT( 13 ) = ULPINV
                  RESULT( 14 ) = ULPINV
                  RESULT( 15 ) = ULPINV
                  GO TO 440
               }

               // Do tests 13 and 14.

               DO 410 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  410          CONTINUE
               DO 420 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  420          CONTINUE
               CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, MAX( 1, M2 ), RESULT( 13 ) )

               NTEST = 15
               DO 430 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  430          CONTINUE
               SRNAMT = 'SSTEVX'
               CALL SSTEVX( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 15 ) = ULPINV
                     GO TO 440
                  }
               }

               // Do test 15.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 15 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

  440          CONTINUE

               NTEST = 16
               DO 450 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  450          CONTINUE
               DO 460 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  460          CONTINUE
               SRNAMT = 'SSTEVD'
               CALL SSTEVD( 'V', N, D1, D2, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVD(V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 16 ) = ULPINV
                     RESULT( 17 ) = ULPINV
                     RESULT( 18 ) = ULPINV
                     GO TO 510
                  }
               }

               // Do tests 16 and 17.

               DO 470 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  470          CONTINUE
               DO 480 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  480          CONTINUE
               CALL SSTT21( N, 0, D3, D4, D1, D2, Z, LDU, WORK, RESULT( 16 ) )

               NTEST = 18
               DO 490 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  490          CONTINUE
               SRNAMT = 'SSTEVD'
               CALL SSTEVD( 'N', N, D3, D4, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVD(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 18 ) = ULPINV
                     GO TO 510
                  }
               }

               // Do test 18.

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 500 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( EVEIGS( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( EVEIGS( J )-D3( J ) ) )
  500          CONTINUE
               RESULT( 18 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  510          CONTINUE

               NTEST = 19
               DO 520 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  520          CONTINUE
               DO 530 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  530          CONTINUE
               SRNAMT = 'SSTEVR'
               CALL SSTEVR( 'V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,I)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 19 ) = ULPINV
                     RESULT( 20 ) = ULPINV
                     RESULT( 21 ) = ULPINV
                     GO TO 570
                  }
               }

               // DO tests 19 and 20.

               DO 540 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  540          CONTINUE
               DO 550 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  550          CONTINUE
               CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, MAX( 1, M2 ), RESULT( 19 ) )


               NTEST = 21
               DO 560 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  560          CONTINUE
               SRNAMT = 'SSTEVR'
               CALL SSTEVR( 'N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,I)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 21 ) = ULPINV
                     GO TO 570
                  }
               }

               // Do test 21.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 21 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )

  570          CONTINUE

               NTEST = 21
               if ( N.GT.0 ) {
                  if ( IL.NE.1 ) {
                     VL = WA1( IL ) - MAX( HALF* ( WA1( IL )-WA1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else {
                     VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
                  if ( IU.NE.N ) {
                     VU = WA1( IU ) + MAX( HALF* ( WA1( IU+1 )-WA1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else {
                     VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
               } else {
                  VL = ZERO
                  VU = ONE
               }

               DO 580 I = 1, N
                  D1( I ) = REAL( A( I, I ) )
  580          CONTINUE
               DO 590 I = 1, N - 1
                  D2( I ) = REAL( A( I+1, I ) )
  590          CONTINUE
               SRNAMT = 'SSTEVR'
               CALL SSTEVR( 'V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 22 ) = ULPINV
                     RESULT( 23 ) = ULPINV
                     RESULT( 24 ) = ULPINV
                     GO TO 630
                  }
               }

               if ( M2.EQ.0 .AND. N.GT.0 ) {
                  RESULT( 22 ) = ULPINV
                  RESULT( 23 ) = ULPINV
                  RESULT( 24 ) = ULPINV
                  GO TO 630
               }

               // Do tests 22 and 23.

               DO 600 I = 1, N
                  D3( I ) = REAL( A( I, I ) )
  600          CONTINUE
               DO 610 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  610          CONTINUE
               CALL SSTT22( N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, MAX( 1, M2 ), RESULT( 22 ) )

               NTEST = 24
               DO 620 I = 1, N - 1
                  D4( I ) = REAL( A( I+1, I ) )
  620          CONTINUE
               SRNAMT = 'SSTEVR'
               CALL SSTEVR( 'N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 24 ) = ULPINV
                     GO TO 630
                  }
               }

               // Do test 24.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( 24 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

  630          CONTINUE



            } else {

               DO 640 I = 1, 24
                  RESULT( I ) = ZERO
  640          CONTINUE
               NTEST = 24
            }

            // Perform remaining tests storing upper or lower triangular
            // part of matrix.

            DO 1720 IUPLO = 0, 1
               if ( IUPLO.EQ.0 ) {
                  UPLO = 'L'
               } else {
                  UPLO = 'U'
               }

               // 4)      Call SSYEV and SSYEVX.

               CALL SLACPY( ' ', N, N, A, LDA, V, LDU )

               NTEST = NTEST + 1
               SRNAMT = 'SSYEV'
               CALL SSYEV( 'V', UPLO, N, A, LDU, D1, WORK, LWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 660
                  }
               }

               // Do tests 25 and 26 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RESULT( NTEST ) )

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               NTEST = NTEST + 2
               SRNAMT = 'SSYEV_2STAGE'
               CALL SSYEV_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK, LWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEV_2STAGE(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 660
                  }
               }

               // Do test 27 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 650 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  650          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  660          CONTINUE
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               NTEST = NTEST + 1

               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
                  if ( IL.NE.1 ) {
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else if ( N.GT.0 ) {
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
                  if ( IU.NE.N ) {
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else if ( N.GT.0 ) {
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
               } else {
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               }

               SRNAMT = 'SSYEVX'
               CALL SSYEVX( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 680
                  }
               }

               // Do tests 28 and 29 (or +54)

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               CALL SSYT21( 1, UPLO, N, 0, A, LDU, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2
               SRNAMT = 'SSYEVX_2STAGE'
               CALL SSYEVX_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVX_2STAGE(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 680
                  }
               }

               // Do test 30 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 670 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  670          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  680          CONTINUE

               NTEST = NTEST + 1
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVX'
               CALL SSYEVX( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 690
                  }
               }

               // Do tests 31 and 32 (or +54)

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVX_2STAGE'
               CALL SSYEVX_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVX_2STAGE(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 690
                  }
               }

               // Do test 33 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
  690          CONTINUE

               NTEST = NTEST + 1
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVX'
               CALL SSYEVX( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 700
                  }
               }

               // Do tests 34 and 35 (or +54)

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVX_2STAGE'
               CALL SSYEVX_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVX_2STAGE(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 700
                  }
               }

               if ( M3.EQ.0 .AND. N.GT.0 ) {
                  RESULT( NTEST ) = ULPINV
                  GO TO 700
               }

               // Do test 36 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

  700          CONTINUE

               // 5)      Call SSPEV and SSPEVX.

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 720 J = 1, N
                     DO 710 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  710                CONTINUE
  720             CONTINUE
               } else {
                  INDX = 1
                  DO 740 J = 1, N
                     DO 730 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  730                CONTINUE
  740             CONTINUE
               }

               NTEST = NTEST + 1
               SRNAMT = 'SSPEV'
               CALL SSPEV( 'V', UPLO, N, WORK, D1, Z, LDU, V, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 800
                  }
               }

               // Do tests 37 and 38 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 760 J = 1, N
                     DO 750 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  750                CONTINUE
  760             CONTINUE
               } else {
                  INDX = 1
                  DO 780 J = 1, N
                     DO 770 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  770                CONTINUE
  780             CONTINUE
               }

               NTEST = NTEST + 2
               SRNAMT = 'SSPEV'
               CALL SSPEV( 'N', UPLO, N, WORK, D3, Z, LDU, V, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 800
                  }
               }

               // Do test 39 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 790 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
  790          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

               // Load array WORK with the upper or lower triangular part
               // of the matrix in packed form.

  800          CONTINUE
               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 820 J = 1, N
                     DO 810 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  810                CONTINUE
  820             CONTINUE
               } else {
                  INDX = 1
                  DO 840 J = 1, N
                     DO 830 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  830                CONTINUE
  840             CONTINUE
               }

               NTEST = NTEST + 1

               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) )
                  if ( IL.NE.1 ) {
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else if ( N.GT.0 ) {
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
                  if ( IU.NE.N ) {
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  } else if ( N.GT.0 ) {
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL )
                  }
               } else {
                  TEMP3 = ZERO
                  VL = ZERO
                  VU = ONE
               }

               SRNAMT = 'SSPEVX'
               CALL SSPEVX( 'V', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 900
                  }
               }

               // Do tests 40 and 41 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 860 J = 1, N
                     DO 850 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  850                CONTINUE
  860             CONTINUE
               } else {
                  INDX = 1
                  DO 880 J = 1, N
                     DO 870 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  870                CONTINUE
  880             CONTINUE
               }

               SRNAMT = 'SSPEVX'
               CALL SSPEVX( 'N', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 900
                  }
               }

               // Do test 42 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 890 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
  890          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

  900          CONTINUE
               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 920 J = 1, N
                     DO 910 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  910                CONTINUE
  920             CONTINUE
               } else {
                  INDX = 1
                  DO 940 J = 1, N
                     DO 930 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  930                CONTINUE
  940             CONTINUE
               }

               NTEST = NTEST + 1

               SRNAMT = 'SSPEVX'
               CALL SSPEVX( 'V', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 990
                  }
               }

               // Do tests 43 and 44 (or +54)

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 960 J = 1, N
                     DO 950 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  950                CONTINUE
  960             CONTINUE
               } else {
                  INDX = 1
                  DO 980 J = 1, N
                     DO 970 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
  970                CONTINUE
  980             CONTINUE
               }

               SRNAMT = 'SSPEVX'
               CALL SSPEVX( 'N', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 990
                  }
               }

               if ( M3.EQ.0 .AND. N.GT.0 ) {
                  RESULT( NTEST ) = ULPINV
                  GO TO 990
               }

               // Do test 45 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

  990          CONTINUE
               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 1010 J = 1, N
                     DO 1000 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1000                CONTINUE
 1010             CONTINUE
               } else {
                  INDX = 1
                  DO 1030 J = 1, N
                     DO 1020 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1020                CONTINUE
 1030             CONTINUE
               }

               NTEST = NTEST + 1

               SRNAMT = 'SSPEVX'
               CALL SSPEVX( 'V', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1080
                  }
               }

               // Do tests 46 and 47 (or +54)

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 1050 J = 1, N
                     DO 1040 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1040                CONTINUE
 1050             CONTINUE
               } else {
                  INDX = 1
                  DO 1070 J = 1, N
                     DO 1060 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1060                CONTINUE
 1070             CONTINUE
               }

               SRNAMT = 'SSPEVX'
               CALL SSPEVX( 'N', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1080
                  }
               }

               if ( M3.EQ.0 .AND. N.GT.0 ) {
                  RESULT( NTEST ) = ULPINV
                  GO TO 1080
               }

               // Do test 48 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

 1080          CONTINUE

               // 6)      Call SSBEV and SSBEVX.

               if ( JTYPE.LE.7 ) {
                  KD = 1
               } else if ( JTYPE.GE.8 .AND. JTYPE.LE.15 ) {
                  KD = MAX( N-1, 0 )
               } else {
                  KD = IHBW
               }

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               if ( IUPLO.EQ.1 ) {
                  DO 1100 J = 1, N
                     DO 1090 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1090                CONTINUE
 1100             CONTINUE
               } else {
                  DO 1120 J = 1, N
                     DO 1110 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1110                CONTINUE
 1120             CONTINUE
               }

               NTEST = NTEST + 1
               SRNAMT = 'SSBEV'
               CALL SSBEV( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1180
                  }
               }

               // Do tests 49 and 50 (or ... )

               CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               if ( IUPLO.EQ.1 ) {
                  DO 1140 J = 1, N
                     DO 1130 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1130                CONTINUE
 1140             CONTINUE
               } else {
                  DO 1160 J = 1, N
                     DO 1150 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1150                CONTINUE
 1160             CONTINUE
               }

               NTEST = NTEST + 2
               SRNAMT = 'SSBEV_2STAGE'
               CALL SSBEV_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSBEV_2STAGE(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1180
                  }
               }

               // Do test 51 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1170 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1170          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

 1180          CONTINUE
               if ( IUPLO.EQ.1 ) {
                  DO 1200 J = 1, N
                     DO 1190 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1190                CONTINUE
 1200             CONTINUE
               } else {
                  DO 1220 J = 1, N
                     DO 1210 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1210                CONTINUE
 1220             CONTINUE
               }

               NTEST = NTEST + 1
               SRNAMT = 'SSBEVX'
               CALL SSBEVX( 'V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1280
                  }
               }

               // Do tests 52 and 53 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2

               if ( IUPLO.EQ.1 ) {
                  DO 1240 J = 1, N
                     DO 1230 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1230                CONTINUE
 1240             CONTINUE
               } else {
                  DO 1260 J = 1, N
                     DO 1250 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1250                CONTINUE
 1260             CONTINUE
               }

               SRNAMT = 'SSBEVX_2STAGE'
               CALL SSBEVX_2STAGE( 'N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSBEVX_2STAGE(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1280
                  }
               }

               // Do test 54 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1270 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA2( J ) ), ABS( WA3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA2( J )-WA3( J ) ) )
 1270          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

 1280          CONTINUE
               NTEST = NTEST + 1
               if ( IUPLO.EQ.1 ) {
                  DO 1300 J = 1, N
                     DO 1290 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1290                CONTINUE
 1300             CONTINUE
               } else {
                  DO 1320 J = 1, N
                     DO 1310 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1310                CONTINUE
 1320             CONTINUE
               }

               SRNAMT = 'SSBEVX'
               CALL SSBEVX( 'V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1370
                  }
               }

               // Do tests 55 and 56 (or +54)

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2

               if ( IUPLO.EQ.1 ) {
                  DO 1340 J = 1, N
                     DO 1330 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1330                CONTINUE
 1340             CONTINUE
               } else {
                  DO 1360 J = 1, N
                     DO 1350 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1350                CONTINUE
 1360             CONTINUE
               }

               SRNAMT = 'SSBEVX_2STAGE'
               CALL SSBEVX_2STAGE( 'N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSBEVX_2STAGE(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1370
                  }
               }

               // Do test 57 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

 1370          CONTINUE
               NTEST = NTEST + 1
               if ( IUPLO.EQ.1 ) {
                  DO 1390 J = 1, N
                     DO 1380 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1380                CONTINUE
 1390             CONTINUE
               } else {
                  DO 1410 J = 1, N
                     DO 1400 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1400                CONTINUE
 1410             CONTINUE
               }

               SRNAMT = 'SSBEVX'
               CALL SSBEVX( 'V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1460
                  }
               }

               // Do tests 58 and 59 (or +54)

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2

               if ( IUPLO.EQ.1 ) {
                  DO 1430 J = 1, N
                     DO 1420 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1420                CONTINUE
 1430             CONTINUE
               } else {
                  DO 1450 J = 1, N
                     DO 1440 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1440                CONTINUE
 1450             CONTINUE
               }

               SRNAMT = 'SSBEVX_2STAGE'
               CALL SSBEVX_2STAGE( 'N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSBEVX_2STAGE(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1460
                  }
               }

               if ( M3.EQ.0 .AND. N.GT.0 ) {
                  RESULT( NTEST ) = ULPINV
                  GO TO 1460
               }

               // Do test 60 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

 1460          CONTINUE

               // 7)      Call SSYEVD

               CALL SLACPY( ' ', N, N, A, LDA, V, LDU )

               NTEST = NTEST + 1
               SRNAMT = 'SSYEVD'
               CALL SSYEVD( 'V', UPLO, N, A, LDU, D1, WORK, LWEDC, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1480
                  }
               }

               // Do tests 61 and 62 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RESULT( NTEST ) )

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               NTEST = NTEST + 2
               SRNAMT = 'SSYEVD_2STAGE'
               CALL SSYEVD_2STAGE( 'N', UPLO, N, A, LDU, D3, WORK,  LWORK, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVD_2STAGE(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1480
                  }
               }

               // Do test 63 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1470 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1470          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

 1480          CONTINUE

               // 8)      Call SSPEVD.

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 1500 J = 1, N
                     DO 1490 I = 1, J
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1490                CONTINUE
 1500             CONTINUE
               } else {
                  INDX = 1
                  DO 1520 J = 1, N
                     DO 1510 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1510                CONTINUE
 1520             CONTINUE
               }

               NTEST = NTEST + 1
               SRNAMT = 'SSPEVD'
               CALL SSPEVD( 'V', UPLO, N, WORK, D1, Z, LDU, WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1580
                  }
               }

               // Do tests 64 and 65 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               if ( IUPLO.EQ.1 ) {
                  INDX = 1
                  DO 1540 J = 1, N
                     DO 1530 I = 1, J

                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1530                CONTINUE
 1540             CONTINUE
               } else {
                  INDX = 1
                  DO 1560 J = 1, N
                     DO 1550 I = J, N
                        WORK( INDX ) = A( I, J )
                        INDX = INDX + 1
 1550                CONTINUE
 1560             CONTINUE
               }

               NTEST = NTEST + 2
               SRNAMT = 'SSPEVD'
               CALL SSPEVD( 'N', UPLO, N, WORK, D3, Z, LDU, WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1580
                  }
               }

               // Do test 66 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1570 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1570          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
 1580          CONTINUE

               // 9)      Call SSBEVD.

               if ( JTYPE.LE.7 ) {
                  KD = 1
               } else if ( JTYPE.GE.8 .AND. JTYPE.LE.15 ) {
                  KD = MAX( N-1, 0 )
               } else {
                  KD = IHBW
               }

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               if ( IUPLO.EQ.1 ) {
                  DO 1600 J = 1, N
                     DO 1590 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1590                CONTINUE
 1600             CONTINUE
               } else {
                  DO 1620 J = 1, N
                     DO 1610 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1610                CONTINUE
 1620             CONTINUE
               }

               NTEST = NTEST + 1
               SRNAMT = 'SSBEVD'
               CALL SSBEVD( 'V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1680
                  }
               }

               // Do tests 67 and 68 (or +54)

               CALL SSYT21( 1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               if ( IUPLO.EQ.1 ) {
                  DO 1640 J = 1, N
                     DO 1630 I = MAX( 1, J-KD ), J
                        V( KD+1+I-J, J ) = A( I, J )
 1630                CONTINUE
 1640             CONTINUE
               } else {
                  DO 1660 J = 1, N
                     DO 1650 I = J, MIN( N, J+KD )
                        V( 1+I-J, J ) = A( I, J )
 1650                CONTINUE
 1660             CONTINUE
               }

               NTEST = NTEST + 2
               SRNAMT = 'SSBEVD_2STAGE'
               CALL SSBEVD_2STAGE( 'N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWORK, IWORK, LIWEDC, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSBEVD_2STAGE(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1680
                  }
               }

               // Do test 69 (or +54)

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1670 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) )
 1670          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

 1680          CONTINUE


               CALL SLACPY( ' ', N, N, A, LDA, V, LDU )
               NTEST = NTEST + 1
               SRNAMT = 'SSYEVR'
               CALL SSYEVR( 'V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1700
                  }
               }

               // Do tests 70 and 71 (or ... )

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               CALL SSYT21( 1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2
               SRNAMT = 'SSYEVR_2STAGE'
               CALL SSYEVR_2STAGE( 'N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVR_2STAGE(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1700
                  }
               }

               // Do test 72 (or ... )

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 1690 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) )
 1690          CONTINUE
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

 1700          CONTINUE

               NTEST = NTEST + 1
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVR'
               CALL SSYEVR( 'V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 1710
                  }
               }

               // Do tests 73 and 74 (or +54)

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVR_2STAGE'
               CALL SSYEVR_2STAGE( 'N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVR_2STAGE(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 1710
                  }
               }

               // Do test 75 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 )
 1710          CONTINUE

               NTEST = NTEST + 1
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVR'
               CALL SSYEVR( 'V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     RESULT( NTEST+1 ) = ULPINV
                     RESULT( NTEST+2 ) = ULPINV
                     GO TO 700
                  }
               }

               // Do tests 76 and 77 (or +54)

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

               CALL SSYT22( 1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) )

               NTEST = NTEST + 2
               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )
               SRNAMT = 'SSYEVR_2STAGE'
               CALL SSYEVR_2STAGE( 'N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 ) 'SSYEVR_2STAGE(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( NTEST ) = ULPINV
                     GO TO 700
                  }
               }

               if ( M3.EQ.0 .AND. N.GT.0 ) {
                  RESULT( NTEST ) = ULPINV
                  GO TO 700
               }

               // Do test 78 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
               if ( N.GT.0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) )
               } else {
                  TEMP3 = ZERO
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

               CALL SLACPY( ' ', N, N, V, LDU, A, LDA )

 1720       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST

            CALL SLAFTS( 'SST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS )

 1730    CONTINUE
 1740 CONTINUE

      // Summary

      CALL ALASVM( 'SST', NOUNIT, NERRS, NTESTT, 0 )

 9999 FORMAT( ' SDRVST2STG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of SDRVST2STG

      }
