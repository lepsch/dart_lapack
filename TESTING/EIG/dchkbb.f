      SUBROUTINE DCHKBB( NSIZES, MVAL, NVAL, NWDTHS, KK, NTYPES, DOTYPE, NRHS, ISEED, THRESH, NOUNIT, A, LDA, AB, LDAB, BD, BE, Q, LDQ, P, LDP, C, LDC, CC, WORK, LWORK, RESULT, INFO )

*  -- LAPACK test routine (input) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDAB, LDC, LDP, LDQ, LWORK, NOUNIT, NRHS, NSIZES, NTYPES, NWDTHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), MVAL( * ), NVAL( * );
      double             A( LDA, * ), AB( LDAB, * ), BD( * ), BE( * ), C( LDC, * ), CC( LDC, * ), P( LDP, * ), Q( LDQ, * ), RESULT( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      int                MAXTYP;
      const              MAXTYP = 15 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KL, KMAX, KU, M, MMAX, MNMAX, MNMIN, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDT01, DBDT02, DGBBRD, DLACPY, DLAHD2, DLASET, DLASUM, DLATMR, DLATMS, DORT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*6, 3*9 /
      DATA               KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3 /
      DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      NTESTT = 0
      INFO = 0

      // Important constants

      BADMM = .FALSE.
      BADNN = .FALSE.
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      DO 10 J = 1, NSIZES
         MMAX = MAX( MMAX, MVAL( J ) )
         IF( MVAL( J ).LT.0 ) BADMM = .TRUE.
         NMAX = MAX( NMAX, NVAL( J ) )
         IF( NVAL( J ).LT.0 ) BADNN = .TRUE.
         MNMAX = MAX( MNMAX, MIN( MVAL( J ), NVAL( J ) ) )
   10 CONTINUE

      BADNNB = .FALSE.
      KMAX = 0
      DO 20 J = 1, NWDTHS
         KMAX = MAX( KMAX, KK( J ) )
         IF( KK( J ).LT.0 ) BADNNB = .TRUE.
   20 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADMM ) {
         INFO = -2
      } else if ( BADNN ) {
         INFO = -3
      } else if ( NWDTHS.LT.0 ) {
         INFO = -4
      } else if ( BADNNB ) {
         INFO = -5
      } else if ( NTYPES.LT.0 ) {
         INFO = -6
      } else if ( NRHS.LT.0 ) {
         INFO = -8
      } else if ( LDA.LT.NMAX ) {
         INFO = -13
      } else if ( LDAB.LT.2*KMAX+1 ) {
         INFO = -15
      } else if ( LDQ.LT.NMAX ) {
         INFO = -19
      } else if ( LDP.LT.NMAX ) {
         INFO = -21
      } else if ( LDC.LT.NMAX ) {
         INFO = -23
      } else if ( ( MAX( LDA, NMAX )+1 )*NMAX.GT.LWORK ) {
         INFO = -26
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DCHKBB', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 .OR. NWDTHS.EQ.0 ) RETURN

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      // Loop over sizes, widths, types

      NERRS = 0
      NMATS = 0

      DO 160 JSIZE = 1, NSIZES
         M = MVAL( JSIZE )
         N = NVAL( JSIZE )
         MNMIN = MIN( M, N )
         AMNINV = ONE / DBLE( MAX( 1, M, N ) )

         DO 150 JWIDTH = 1, NWDTHS
            K = KK( JWIDTH )
            IF( K.GE.M .AND. K.GE.N ) GO TO 150
            KL = MAX( 0, MIN( M-1, K ) )
            KU = MAX( 0, MIN( N-1, K ) )

            if ( NSIZES.NE.1 ) {
               MTYPES = MIN( MAXTYP, NTYPES )
            } else {
               MTYPES = MIN( MAXTYP+1, NTYPES )
            }

            DO 140 JTYPE = 1, MTYPES
               IF( .NOT.DOTYPE( JTYPE ) ) GO TO 140
               NMATS = NMATS + 1
               NTEST = 0

               DO 30 J = 1, 4
                  IOLDSD( J ) = ISEED( J )
   30          CONTINUE

               // Compute "A".

               // Control parameters:

                   // KMAGN  KMODE        KTYPE
               // =1  O(1)   clustered 1  zero
               // =2  large  clustered 2  identity
               // =3  small  exponential  (none)
               // =4         arithmetic   diagonal, (w/ singular values)
               // =5         random log   (none)
               // =6         random       nonhermitian, w/ singular values
               // =7                      (none)
               // =8                      (none)
               // =9                      random nonhermitian

               IF( MTYPES.GT.MAXTYP ) GO TO 90

               ITYPE = KTYPE( JTYPE )
               IMODE = KMODE( JTYPE )

               // Compute norm

               GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40          CONTINUE
               ANORM = ONE
               GO TO 70

   50          CONTINUE
               ANORM = ( RTOVFL*ULP )*AMNINV
               GO TO 70

   60          CONTINUE
               ANORM = RTUNFL*MAX( M, N )*ULPINV
               GO TO 70

   70          CONTINUE

               CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               CALL DLASET( 'Full', LDAB, N, ZERO, ZERO, AB, LDAB )
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
   80             CONTINUE

               } else if ( ITYPE.EQ.4 ) {

                  // Diagonal Matrix, singular values specified

                  CALL DLATMS( M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( M+1 ), IINFO )

               } else if ( ITYPE.EQ.6 ) {

                  // Nonhermitian, singular values specified

                  CALL DLATMS( M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, KL, KU, 'N', A, LDA, WORK( M+1 ), IINFO )

               } else if ( ITYPE.EQ.9 ) {

                  // Nonhermitian, random entries

                  CALL DLATMR( M, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, KL, KU, ZERO, ANORM, 'N', A, LDA, IDUMMA, IINFO )

               } else {

                  IINFO = 1
               }

               // Generate Right-Hand Side

               CALL DLATMR( M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( M+1 ), 1, ONE, WORK( 2*M+1 ), 1, ONE, 'N', IDUMMA, M, NRHS, ZERO, ONE, 'NO', C, LDC, IDUMMA, IINFO )

               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

   90          CONTINUE

               // Copy A to band storage.

               DO 110 J = 1, N
                  DO 100 I = MAX( 1, J-KU ), MIN( M, J+KL )
                     AB( KU+1+I-J, J ) = A( I, J )
  100             CONTINUE
  110          CONTINUE

               // Copy C

               CALL DLACPY( 'Full', M, NRHS, C, LDC, CC, LDC )

               // Call DGBBRD to compute B, Q and P, and to update C.

               CALL DGBBRD( 'B', M, N, NRHS, KL, KU, AB, LDAB, BD, BE, Q, LDQ, P, LDP, CC, LDC, WORK, IINFO )

               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DGBBRD', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 1 ) = ULPINV
                     GO TO 120
                  }
               }

               // Test 1:  Check the decomposition A := Q * B * P'
                    // 2:  Check the orthogonality of Q
                    // 3:  Check the orthogonality of P
                    // 4:  Check the computation of Q' * C

               CALL DBDT01( M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RESULT( 1 ) )                CALL DORT01( 'Columns', M, M, Q, LDQ, WORK, LWORK, RESULT( 2 ) )                CALL DORT01( 'Rows', N, N, P, LDP, WORK, LWORK, RESULT( 3 ) )                CALL DBDT02( M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RESULT( 4 ) )

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 4
  120          CONTINUE
               NTESTT = NTESTT + NTEST

               // Print out tests which fail.

               DO 130 JR = 1, NTEST
                  if ( RESULT( JR ).GE.THRESH ) {
                     IF( NERRS.EQ.0 ) CALL DLAHD2( NOUNIT, 'DBB' )
                     NERRS = NERRS + 1
                     WRITE( NOUNIT, FMT = 9998 )M, N, K, IOLDSD, JTYPE, JR, RESULT( JR )
                  }
  130          CONTINUE

  140       CONTINUE
  150    CONTINUE
  160 CONTINUE

      // Summary

      CALL DLASUM( 'DBB', NOUNIT, NERRS, NTESTT )
      RETURN

 9999 FORMAT( ' DCHKBB: ', A, ' returned INFO=', I5, '.', / 9X, 'M=',
     $      I5, ' N=', I5, ' K=', I5, ', JTYPE=', I5, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' M =', I4, ' N=', I4, ', K=', I3, ', seed=',
     $      4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 )

      // End of DCHKBB

      }
