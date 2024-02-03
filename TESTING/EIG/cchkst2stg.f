      SUBROUTINE CCHKST2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, AP, SD, SE, D1, D2, D3, D4, D5, WA1, WA2, WA3, WR, U, LDU, V, VP, TAU, Z, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LRWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               D1( * ), D2( * ), D3( * ), D4( * ), D5( * ), RESULT( * ), RWORK( * ), SD( * ), SE( * ), WA1( * ), WA2( * ), WA3( * ), WR( * )
      COMPLEX            A( LDA, * ), AP( * ), TAU( * ), U( LDU, * ), V( LDU, * ), VP( * ), WORK( * ), Z( LDU, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, EIGHT, TEN, HUN
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, EIGHT = 8.0E0, TEN = 10.0E0, HUN = 100.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      REAL               HALF
      const              HALF = ONE / TWO ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      bool               CRANGE;
      const              CRANGE = .FALSE. ;
      bool               CREL;
      const              CREL = .FALSE. ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, TRYRAC;
      int                I, IINFO, IL, IMODE, INDE, INDRWK, ITEMP, ITYPE, IU, J, JC, JR, JSIZE, JTYPE, LGN, LIWEDC, LOG2UI, LRWEDC, LWEDC, M, M2, M3, MTYPES, N, NAP, NBLOCK, NERRS, NMATS, NMAX, NSPLIT, NTEST, NTESTT, LH, LW;
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, TEMP4, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      REAL               DUMMA( 1 )
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SLAMCH, SLARND, SSXT1
      // EXTERNAL ILAENV, SLAMCH, SLARND, SSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLASUM, SSTEBZ, SSTECH, SSTERF, XERBLA, CCOPY, CHET21, CHETRD, CHPT21, CHPTRD, CLACPY, CLASET, CLATMR, CLATMS, CPTEQR, CSTEDC, CSTEMR, CSTEIN, CSTEQR, CSTT21, CSTT22, CUNGTR, CUPGTR, CHETRD_2STAGE, SLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, CONJG, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 10 /       DATA               KMAGN / 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1 /       DATA               KMODE / 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 3, 1, 4, 4, 3 /
      // ..
      // .. Executable Statements ..

      // Keep ftnchek happy
      IDUMMA( 1 ) = 1

      // Check for errors

      NTESTT = 0
      INFO = 0

      // Important constants

      BADNN = .FALSE.
      TRYRAC = .TRUE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      NBLOCK = ILAENV( 1, 'CHETRD', 'L', NMAX, -1, -1, -1 )
      NBLOCK = MIN( NMAX, MAX( 1, NBLOCK ) )

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
         INFO = -23
      } else if ( 2*MAX( 2, NMAX )**2.GT.LWORK ) {
         INFO = -29
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CCHKST2STG', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) )
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      // Loop over sizes, types

      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
   20 CONTINUE
      NERRS = 0
      NMATS = 0

      DO 310 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         if ( N.GT.0 ) {
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2
            LRWEDC = 1 + 3*N + 2*N*LGN + 4*N**2
            LIWEDC = 6 + 6*N + 5*N*LGN
         } else {
            LWEDC = 8
            LRWEDC = 7
            LIWEDC = 12
         }
         NAP = ( N*( N+1 ) ) / 2
         ANINV = ONE / REAL( MAX( 1, N ) )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 300 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 300
            NMATS = NMATS + 1
            NTEST = 0

            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE

            // Compute "A"

            // Control parameters:

                // KMAGN  KMODE        KTYPE
            // =1  O(1)   clustered 1  zero
            // =2  large  clustered 2  identity
            // =3  small  exponential  (none)
            // =4         arithmetic   diagonal, (w/ eigenvalues)
            // =5         random log   Hermitian, w/ eigenvalues
            // =6         random       (none)
            // =7                      random diagonal
            // =8                      random Hermitian
            // =9                      positive definite
            // =10                     diagonally dominant tridiagonal

            IF( MTYPES.GT.MAXTYP ) GO TO 100

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

            CALL CLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
            IINFO = 0
            if ( JTYPE.LE.15 ) {
               COND = ULPINV
            } else {
               COND = ULPINV*ANINV / TEN
            }

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE.EQ.1 ) {
               IINFO = 0

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               DO 80 JC = 1, N
                  A( JC, JC ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO )


            } else if ( ITYPE.EQ.5 ) {

               // Hermitian, eigenvalues specified

               CALL CLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO )

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               CALL CLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.8 ) {

               // Hermitian, random eigenvalues

               CALL CLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.9 ) {

               // Positive definite, eigenvalues specified.

               CALL CLATMS( N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO )

            } else if ( ITYPE.EQ.10 ) {

               // Positive definite tridiagonal, eigenvalues specified.

               CALL CLATMS( N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, 1, 1, 'N', A, LDA, WORK, IINFO )
               DO 90 I = 2, N
                  TEMP1 = ABS( A( I-1, I ) )
                  TEMP2 = SQRT( ABS( A( I-1, I-1 )*A( I, I ) ) )
                  if ( TEMP1.GT.HALF*TEMP2 ) {
                     A( I-1, I ) = A( I-1, I )* ( HALF*TEMP2 / ( UNFL+TEMP1 ) )
                     A( I, I-1 ) = CONJG( A( I-1, I ) )
                  }
   90          CONTINUE

            } else {

               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

  100       CONTINUE

            // Call CHETRD and CUNGTR to compute S and U from
            // upper triangle.

            CALL CLACPY( 'U', N, N, A, LDA, V, LDU )

            NTEST = 1
            CALL CHETRD( 'U', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO )

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHETRD(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 1 ) = ULPINV
                  GO TO 280
               }
            }

            CALL CLACPY( 'U', N, N, V, LDU, U, LDU )

            NTEST = 2
            CALL CUNGTR( 'U', N, U, LDU, TAU, WORK, LWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUNGTR(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 2 ) = ULPINV
                  GO TO 280
               }
            }

            // Do tests 1 and 2

            CALL CHET21( 2, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RWORK, RESULT( 1 ) )             CALL CHET21( 3, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RWORK, RESULT( 2 ) )

            // Compute D1 the eigenvalues resulting from the tridiagonal
            // form using the standard 1-stage algorithm and use it as a
            // reference to compare with the 2-stage technique

            // Compute D1 from the 1-stage and used as reference for the
            // 2-stage

            CALL SCOPY( N, SD, 1, D1, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

            CALL CSTEQR( 'N', N, D1, RWORK, WORK, LDU, RWORK( N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 3 ) = ULPINV
                  GO TO 280
               }
            }

            // 2-STAGE TRD Upper case is used to compute D2.
            // Note to set SD and SE to zero to be sure not reusing
            // the one from above. Compare it with D1 computed
            // using the 1-stage.

            CALL SLASET( 'Full', N, 1, ZERO, ZERO, SD, N )
            CALL SLASET( 'Full', N, 1, ZERO, ZERO, SE, N )
            CALL CLACPY( 'U', N, N, A, LDA, V, LDU )
            LH = MAX(1, 4*N)
            LW = LWORK - LH
            CALL CHETRD_2STAGE( 'N', "U", N, V, LDU, SD, SE, TAU,  WORK, LH, WORK( LH+1 ), LW, IINFO )

            // Compute D2 from the 2-stage Upper case

            CALL SCOPY( N, SD, 1, D2, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

            NTEST = 3
            CALL CSTEQR( 'N', N, D2, RWORK, WORK, LDU, RWORK( N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 3 ) = ULPINV
                  GO TO 280
               }
            }

            // 2-STAGE TRD Lower case is used to compute D3.
            // Note to set SD and SE to zero to be sure not reusing
            // the one from above. Compare it with D1 computed
            // using the 1-stage.

            CALL SLASET( 'Full', N, 1, ZERO, ZERO, SD, N )
            CALL SLASET( 'Full', N, 1, ZERO, ZERO, SE, N )
            CALL CLACPY( 'L', N, N, A, LDA, V, LDU )
            CALL CHETRD_2STAGE( 'N', "L", N, V, LDU, SD, SE, TAU,  WORK, LH, WORK( LH+1 ), LW, IINFO )

            // Compute D3 from the 2-stage Upper case

            CALL SCOPY( N, SD, 1, D3, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

            NTEST = 4
            CALL CSTEQR( 'N', N, D3, RWORK, WORK, LDU, RWORK( N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 4 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 3 and 4 which are similar to 11 and 12 but with the
            // D1 computed using the standard 1-stage reduction as reference

            NTEST = 4
            TEMP1 = ZERO
            TEMP2 = ZERO
            TEMP3 = ZERO
            TEMP4 = ZERO

            DO 151 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
               TEMP3 = MAX( TEMP3, ABS( D1( J ) ), ABS( D3( J ) ) )
               TEMP4 = MAX( TEMP4, ABS( D1( J )-D3( J ) ) )
  151       CONTINUE

            RESULT( 3 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
            RESULT( 4 ) = TEMP4 / MAX( UNFL, ULP*MAX( TEMP3, TEMP4 ) )

            // Store the upper triangle of A in AP

            I = 0
            DO 120 JC = 1, N
               DO 110 JR = 1, JC
                  I = I + 1
                  AP( I ) = A( JR, JC )
  110          CONTINUE
  120       CONTINUE

            // Call CHPTRD and CUPGTR to compute S and U from AP

            CALL CCOPY( NAP, AP, 1, VP, 1 )

            NTEST = 5
            CALL CHPTRD( 'U', N, VP, SD, SE, TAU, IINFO )

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHPTRD(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 5 ) = ULPINV
                  GO TO 280
               }
            }

            NTEST = 6
            CALL CUPGTR( 'U', N, VP, TAU, U, LDU, WORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUPGTR(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 6 ) = ULPINV
                  GO TO 280
               }
            }

            // Do tests 5 and 6

            CALL CHPT21( 2, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 5 ) )             CALL CHPT21( 3, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 6 ) )

            // Store the lower triangle of A in AP

            I = 0
            DO 140 JC = 1, N
               DO 130 JR = JC, N
                  I = I + 1
                  AP( I ) = A( JR, JC )
  130          CONTINUE
  140       CONTINUE

            // Call CHPTRD and CUPGTR to compute S and U from AP

            CALL CCOPY( NAP, AP, 1, VP, 1 )

            NTEST = 7
            CALL CHPTRD( 'L', N, VP, SD, SE, TAU, IINFO )

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHPTRD(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 7 ) = ULPINV
                  GO TO 280
               }
            }

            NTEST = 8
            CALL CUPGTR( 'L', N, VP, TAU, U, LDU, WORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUPGTR(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 8 ) = ULPINV
                  GO TO 280
               }
            }

            CALL CHPT21( 2, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 7 ) )             CALL CHPT21( 3, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 8 ) )

            // Call CSTEQR to compute D1, D2, and Z, do tests.

            // Compute D1 and Z

            CALL SCOPY( N, SD, 1, D1, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )
            CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

            NTEST = 9
            CALL CSTEQR( 'V', N, D1, RWORK, Z, LDU, RWORK( N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 9 ) = ULPINV
                  GO TO 280
               }
            }

            // Compute D2

            CALL SCOPY( N, SD, 1, D2, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

            NTEST = 11
            CALL CSTEQR( 'N', N, D2, RWORK, WORK, LDU, RWORK( N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 11 ) = ULPINV
                  GO TO 280
               }
            }

            // Compute D3 (using PWK method)

            CALL SCOPY( N, SD, 1, D3, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

            NTEST = 12
            CALL SSTERF( N, D3, RWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTERF', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 12 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 9 and 10

            CALL CSTT21( N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 9 ) )

            // Do Tests 11 and 12

            TEMP1 = ZERO
            TEMP2 = ZERO
            TEMP3 = ZERO
            TEMP4 = ZERO

            DO 150 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
               TEMP3 = MAX( TEMP3, ABS( D1( J ) ), ABS( D3( J ) ) )
               TEMP4 = MAX( TEMP4, ABS( D1( J )-D3( J ) ) )
  150       CONTINUE

            RESULT( 11 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
            RESULT( 12 ) = TEMP4 / MAX( UNFL, ULP*MAX( TEMP3, TEMP4 ) )

            // Do Test 13 -- Sturm Sequence Test of Eigenvalues
                          // Go up by factors of two until it succeeds

            NTEST = 13
            TEMP1 = THRESH*( HALF-ULP )

            DO 160 J = 0, LOG2UI
               CALL SSTECH( N, SD, SE, D1, TEMP1, RWORK, IINFO )
               IF( IINFO.EQ.0 ) GO TO 170
               TEMP1 = TEMP1*TWO
  160       CONTINUE

  170       CONTINUE
            RESULT( 13 ) = TEMP1

            // For positive definite matrices ( JTYPE.GT.15 ) call CPTEQR
            // and do tests 14, 15, and 16 .

            if ( JTYPE.GT.15 ) {

               // Compute D4 and Z4

               CALL SCOPY( N, SD, 1, D4, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )
               CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

               NTEST = 14
               CALL CPTEQR( 'V', N, D4, RWORK, Z, LDU, RWORK( N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CPTEQR(V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 14 ) = ULPINV
                     GO TO 280
                  }
               }

               // Do Tests 14 and 15

               CALL CSTT21( N, 0, SD, SE, D4, DUMMA, Z, LDU, WORK, RWORK, RESULT( 14 ) )

               // Compute D5

               CALL SCOPY( N, SD, 1, D5, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

               NTEST = 16
               CALL CPTEQR( 'N', N, D5, RWORK, Z, LDU, RWORK( N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CPTEQR(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 16 ) = ULPINV
                     GO TO 280
                  }
               }

               // Do Test 16

               TEMP1 = ZERO
               TEMP2 = ZERO
               DO 180 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D4( J ) ), ABS( D5( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D4( J )-D5( J ) ) )
  180          CONTINUE

               RESULT( 16 ) = TEMP2 / MAX( UNFL, HUN*ULP*MAX( TEMP1, TEMP2 ) )
            } else {
               RESULT( 14 ) = ZERO
               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
            }

            // Call SSTEBZ with different options and do tests 17-18.

               // If S is positive definite and diagonally dominant,
               // ask for all eigenvalues with high relative accuracy.

            VL = ZERO
            VU = ZERO
            IL = 0
            IU = 0
            if ( JTYPE.EQ.21 ) {
               NTEST = 17
               ABSTOL = UNFL + UNFL
               CALL SSTEBZ( 'A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WR, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(A,rel)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 17 ) = ULPINV
                     GO TO 280
                  }
               }

               // Do test 17

               TEMP2 = TWO*( TWO*N-ONE )*ULP*( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4

               TEMP1 = ZERO
               DO 190 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D4( J )-WR( N-J+1 ) ) / ( ABSTOL+ABS( D4( J ) ) ) )
  190          CONTINUE

               RESULT( 17 ) = TEMP1 / TEMP2
            } else {
               RESULT( 17 ) = ZERO
            }

            // Now ask for all eigenvalues with high absolute accuracy.

            NTEST = 18
            ABSTOL = UNFL + UNFL
            CALL SSTEBZ( 'A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(A)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 18 ) = ULPINV
                  GO TO 280
               }
            }

            // Do test 18

            TEMP1 = ZERO
            TEMP2 = ZERO
            DO 200 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( D3( J ) ), ABS( WA1( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D3( J )-WA1( J ) ) )
  200       CONTINUE

            RESULT( 18 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

            // Choose random values for IL and IU, and ask for the
            // IL-th through IU-th eigenvalues.

            NTEST = 19
            if ( N.LE.1 ) {
               IL = 1
               IU = N
            } else {
               IL = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
               IU = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
               if ( IU.LT.IL ) {
                  ITEMP = IU
                  IU = IL
                  IL = ITEMP
               }
            }

            CALL SSTEBZ( 'I', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M2, NSPLIT, WA2, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(I)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 19 ) = ULPINV
                  GO TO 280
               }
            }

            // Determine the values VL and VU of the IL-th and IU-th
            // eigenvalues and ask for all eigenvalues in this range.

            if ( N.GT.0 ) {
               if ( IL.NE.1 ) {
                  VL = WA1( IL ) - MAX( HALF*( WA1( IL )-WA1( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
               } else {
                  VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ), ULP*ANORM, TWO*RTUNFL )
               }
               if ( IU.NE.N ) {
                  VU = WA1( IU ) + MAX( HALF*( WA1( IU+1 )-WA1( IU ) ), ULP*ANORM, TWO*RTUNFL )
               } else {
                  VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ), ULP*ANORM, TWO*RTUNFL )
               }
            } else {
               VL = ZERO
               VU = ONE
            }

            CALL SSTEBZ( 'V', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M3, NSPLIT, WA3, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 19 ) = ULPINV
                  GO TO 280
               }
            }

            if ( M3.EQ.0 .AND. N.NE.0 ) {
               RESULT( 19 ) = ULPINV
               GO TO 280
            }

            // Do test 19

            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            if ( N.GT.0 ) {
               TEMP3 = MAX( ABS( WA1( N ) ), ABS( WA1( 1 ) ) )
            } else {
               TEMP3 = ZERO
            }

            RESULT( 19 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

            // Call CSTEIN to compute eigenvectors corresponding to
            // eigenvalues in WA1.  (First call SSTEBZ again, to make sure
            // it returns these eigenvalues in the correct order.)

            NTEST = 21
            CALL SSTEBZ( 'A', 'B', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(A,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 20 ) = ULPINV
                  RESULT( 21 ) = ULPINV
                  GO TO 280
               }
            }

            CALL CSTEIN( N, SD, SE, M, WA1, IWORK( 1 ), IWORK( N+1 ), Z, LDU, RWORK, IWORK( 2*N+1 ), IWORK( 3*N+1 ), IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEIN', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 20 ) = ULPINV
                  RESULT( 21 ) = ULPINV
                  GO TO 280
               }
            }

            // Do tests 20 and 21

            CALL CSTT21( N, 0, SD, SE, WA1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 20 ) )

            // Call CSTEDC(I) to compute D1 and Z, do tests.

            // Compute D1 and Z

            INDE = 1
            INDRWK = INDE + N
            CALL SCOPY( N, SD, 1, D1, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK( INDE ), 1 )
            CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

            NTEST = 22
            CALL CSTEDC( 'I', N, D1, RWORK( INDE ), Z, LDU, WORK, LWEDC, RWORK( INDRWK ), LRWEDC, IWORK, LIWEDC, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEDC(I)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 22 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 22 and 23

            CALL CSTT21( N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 22 ) )

            // Call CSTEDC(V) to compute D1 and Z, do tests.

            // Compute D1 and Z

            CALL SCOPY( N, SD, 1, D1, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK( INDE ), 1 )
            CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

            NTEST = 24
            CALL CSTEDC( 'V', N, D1, RWORK( INDE ), Z, LDU, WORK, LWEDC, RWORK( INDRWK ), LRWEDC, IWORK, LIWEDC, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEDC(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 24 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 24 and 25

            CALL CSTT21( N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 24 ) )

            // Call CSTEDC(N) to compute D2, do tests.

            // Compute D2

            CALL SCOPY( N, SD, 1, D2, 1 )
            IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK( INDE ), 1 )
            CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

            NTEST = 26
            CALL CSTEDC( 'N', N, D2, RWORK( INDE ), Z, LDU, WORK, LWEDC, RWORK( INDRWK ), LRWEDC, IWORK, LIWEDC, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEDC(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 26 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Test 26

            TEMP1 = ZERO
            TEMP2 = ZERO

            DO 210 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
  210       CONTINUE

            RESULT( 26 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

            // Only test CSTEMR if IEEE compliant

            if ( ILAENV( 10, 'CSTEMR', 'VA', 1, 0, 0, 0 ).EQ.1 .AND. ILAENV( 11, 'CSTEMR', 'VA', 1, 0, 0, 0 ).EQ.1 ) {

            // Call CSTEMR, do test 27 (relative eigenvalue accuracy)

               // If S is positive definite and diagonally dominant,
               // ask for all eigenvalues with high relative accuracy.

               VL = ZERO
               VU = ZERO
               IL = 0
               IU = 0
               if ( JTYPE.EQ.21 .AND. CREL ) {
                  NTEST = 27
                  ABSTOL = UNFL + UNFL
                  CALL CSTEMR( 'V', 'A', N, SD, SE, VL, VU, IL, IU, M, WR, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK, LRWORK, IWORK( 2*N+1 ), LWORK-2*N, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,A,rel)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 27 ) = ULPINV
                        GO TO 270
                     }
                  }

               // Do test 27

                  TEMP2 = TWO*( TWO*N-ONE )*ULP*( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4

                  TEMP1 = ZERO
                  DO 220 J = 1, N
                     TEMP1 = MAX( TEMP1, ABS( D4( J )-WR( N-J+1 ) ) / ( ABSTOL+ABS( D4( J ) ) ) )
  220             CONTINUE

                  RESULT( 27 ) = TEMP1 / TEMP2

                  IL = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  IU = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  if ( IU.LT.IL ) {
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  }

                  if ( CRANGE ) {
                     NTEST = 28
                     ABSTOL = UNFL + UNFL
                     CALL CSTEMR( 'V', 'I', N, SD, SE, VL, VU, IL, IU, M, WR, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK, LRWORK, IWORK( 2*N+1 ), LWORK-2*N, IINFO )

                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,I,rel)', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( 28 ) = ULPINV
                           GO TO 270
                        }
                     }

                  // Do test 28

                     TEMP2 = TWO*( TWO*N-ONE )*ULP* ( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4

                     TEMP1 = ZERO
                     DO 230 J = IL, IU
                        TEMP1 = MAX( TEMP1, ABS( WR( J-IL+1 )-D4( N-J+ 1 ) ) / ( ABSTOL+ABS( WR( J-IL+1 ) ) ) )
  230                CONTINUE

                     RESULT( 28 ) = TEMP1 / TEMP2
                  } else {
                     RESULT( 28 ) = ZERO
                  }
               } else {
                  RESULT( 27 ) = ZERO
                  RESULT( 28 ) = ZERO
               }

            // Call CSTEMR(V,I) to compute D1 and Z, do tests.

            // Compute D1 and Z

               CALL SCOPY( N, SD, 1, D5, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )
               CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

               if ( CRANGE ) {
                  NTEST = 29
                  IL = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  IU = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  if ( IU.LT.IL ) {
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  }
                  CALL CSTEMR( 'V', 'I', N, D5, RWORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,I)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 29 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Tests 29 and 30

            // Call CSTEMR to compute D2, do tests.

            // Compute D2

                  CALL SCOPY( N, SD, 1, D5, 1 )
                  IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

                  NTEST = 31
                  CALL CSTEMR( 'N', 'I', N, D5, RWORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(N,I)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 31 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Test 31

                  TEMP1 = ZERO
                  TEMP2 = ZERO

                  DO 240 J = 1, IU - IL + 1
                     TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                     TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
  240             CONTINUE

                  RESULT( 31 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

            // Call CSTEMR(V,V) to compute D1 and Z, do tests.

            // Compute D1 and Z

                  CALL SCOPY( N, SD, 1, D5, 1 )
                  IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )
                  CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDU )

                  NTEST = 32

                  if ( N.GT.0 ) {
                     if ( IL.NE.1 ) {
                        VL = D2( IL ) - MAX( HALF* ( D2( IL )-D2( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
                     } else {
                        VL = D2( 1 ) - MAX( HALF*( D2( N )-D2( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                     }
                     if ( IU.NE.N ) {
                        VU = D2( IU ) + MAX( HALF* ( D2( IU+1 )-D2( IU ) ), ULP*ANORM, TWO*RTUNFL )
                     } else {
                        VU = D2( N ) + MAX( HALF*( D2( N )-D2( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                     }
                  } else {
                     VL = ZERO
                     VU = ONE
                  }

                  CALL CSTEMR( 'V', 'V', N, D5, RWORK, VL, VU, IL, IU, M, D1, Z, LDU, M, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,V)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 32 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Tests 32 and 33

                  CALL CSTT22( N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RWORK, RESULT( 32 ) )

            // Call CSTEMR to compute D2, do tests.

            // Compute D2

                  CALL SCOPY( N, SD, 1, D5, 1 )
                  IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

                  NTEST = 34
                  CALL CSTEMR( 'N', 'V', N, D5, RWORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(N,V)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 34 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Test 34

                  TEMP1 = ZERO
                  TEMP2 = ZERO

                  DO 250 J = 1, IU - IL + 1
                     TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                     TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
  250             CONTINUE

                  RESULT( 34 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
               } else {
                  RESULT( 29 ) = ZERO
                  RESULT( 30 ) = ZERO
                  RESULT( 31 ) = ZERO
                  RESULT( 32 ) = ZERO
                  RESULT( 33 ) = ZERO
                  RESULT( 34 ) = ZERO
               }

            // Call CSTEMR(V,A) to compute D1 and Z, do tests.

            // Compute D1 and Z

               CALL SCOPY( N, SD, 1, D5, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

               NTEST = 35

               CALL CSTEMR( 'V', 'A', N, D5, RWORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 35 ) = ULPINV
                     GO TO 280
                  }
               }

            // Do Tests 35 and 36

               CALL CSTT22( N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RWORK, RESULT( 35 ) )

            // Call CSTEMR to compute D2, do tests.

            // Compute D2

               CALL SCOPY( N, SD, 1, D5, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, RWORK, 1 )

               NTEST = 37
               CALL CSTEMR( 'N', 'A', N, D5, RWORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CSTEMR(N,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 37 ) = ULPINV
                     GO TO 280
                  }
               }

            // Do Test 37

               TEMP1 = ZERO
               TEMP2 = ZERO

               DO 260 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
  260          CONTINUE

               RESULT( 37 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
            }
  270       CONTINUE
  280       CONTINUE
            NTESTT = NTESTT + NTEST

            // End of Loop -- Check for RESULT(j) > THRESH

            // Print out tests which fail.

            DO 290 JR = 1, NTEST
               if ( RESULT( JR ).GE.THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS.EQ.0 ) {
                     WRITE( NOUNIT, FMT = 9998 )'CST'
                     WRITE( NOUNIT, FMT = 9997 )
                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )'Hermitian'
                     WRITE( NOUNIT, FMT = 9994 )

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9987 )
                  }
                  NERRS = NERRS + 1
                  if ( RESULT( JR ).LT.10000.0E0 ) {
                     WRITE( NOUNIT, FMT = 9989 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  } else {
                     WRITE( NOUNIT, FMT = 9988 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  }
               }
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE

      // Summary

      CALL SLASUM( 'CST', NOUNIT, NERRS, NTESTT )
      RETURN

 9999 FORMAT( ' CCHKST2STG: ', A, ' returned INFO=', I6, '.', / 9X,
     $   'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9998 FORMAT( / 1X, A3, ' -- Complex Hermitian eigenvalue problem' )
 9997 FORMAT( ' Matrix types (see CCHKST2STG for details): ' )

 9996 FORMAT( / ' Special Matrices:',
     $      / '  1=Zero matrix.                        ',
     $      '  5=Diagonal: clustered entries.',
     $      / '  2=Identity matrix.                    ',
     $      '  6=Diagonal: large, evenly spaced.',
     $      / '  3=Diagonal: evenly spaced entries.    ',
     $      '  7=Diagonal: small, evenly spaced.',
     $      / '  4=Diagonal: geometr. spaced entries.' )
 9995 FORMAT( ' Dense ', A, ' Matrices:',
     $      / '  8=Evenly spaced eigenvals.            ',
     $      ' 12=Small, evenly spaced eigenvals.',
     $      / '  9=Geometrically spaced eigenvals.     ',
     $      ' 13=Matrix with random O(1) entries.',
     $      / ' 10=Clustered eigenvalues.              ',
     $      ' 14=Matrix with large random entries.',
     $      / ' 11=Large, evenly spaced eigenvals.     ',
     $      ' 15=Matrix with small random entries.' )
 9994 FORMAT( ' 16=Positive definite, evenly spaced eigenvalues',
     $      / ' 17=Positive definite, geometrically spaced eigenvlaues',
     $      / ' 18=Positive definite, clustered eigenvalues',
     $      / ' 19=Positive definite, small evenly spaced eigenvalues',
     $      / ' 20=Positive definite, large evenly spaced eigenvalues',
     $      / ' 21=Diagonally dominant tridiagonal, geometrically',
     $      ' spaced eigenvalues' )

 9989 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I3, ' is', 0P, F8.2 )
 9988 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I3, ' is', 1P, E10.3 )

 9987 FORMAT( / 'Test performed:  see CCHKST2STG for details.', / )

      // End of CCHKST2STG

      }
