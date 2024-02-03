      SUBROUTINE ZDRVRF2( NOUT, NN, NVAL, A, LDA, ARF, AP, ASAV  )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, NN, NOUT;
*     ..
*     .. Array Arguments ..
      int                NVAL( NN );
      COMPLEX*16         A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * )
*     ..
*
*  =====================================================================
*     ..
*     .. Local Scalars ..
      bool               LOWER, OK1, OK2;
      String             UPLO, CFORM;
      int                I, IFORM, IIN, INFO, IUPLO, J, N, NERRS, NRUN;
*     ..
*     .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
*     ..
*     .. External Functions ..
      COMPLEX*16         ZLARND
      // EXTERNAL ZLARND
*     ..
*     .. External Subroutines ..
      // EXTERNAL ZTFTTR, ZTFTTP, ZTRTTF, ZTRTTP, ZTPTTR, ZTPTTF
*     ..
*     .. Scalars in Common ..
      String             SRNAMT;
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FORMS / 'N', 'C' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      NRUN = 0
      NERRS = 0
      INFO = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      DO 120 IIN = 1, NN
*
         N = NVAL( IIN )
*
*        Do first for UPLO = 'U', then for UPLO = 'L'
*
         DO 110 IUPLO = 1, 2
*
            UPLO = UPLOS( IUPLO )
            LOWER = .TRUE.
            IF ( IUPLO.EQ.1 ) LOWER = .FALSE.
*
*           Do first for CFORM = 'N', then for CFORM = 'C'
*
            DO 100 IFORM = 1, 2
*
               CFORM = FORMS( IFORM )
*
               NRUN = NRUN + 1
*
               DO J = 1, N
                  DO I = 1, N
                     A( I, J) = ZLARND( 4, ISEED )
                  END DO
               END DO
*
               SRNAMT = 'ZTRTTF'
               CALL ZTRTTF( CFORM, UPLO, N, A, LDA, ARF, INFO )
*
               SRNAMT = 'ZTFTTP'
               CALL ZTFTTP( CFORM, UPLO, N, ARF, AP, INFO )
*
               SRNAMT = 'ZTPTTR'
               CALL ZTPTTR( UPLO, N, AP, ASAV, LDA, INFO )
*
               OK1 = .TRUE.
               IF ( LOWER ) THEN
                  DO J = 1, N
                     DO I = J, N
                        IF ( A(I,J).NE.ASAV(I,J) ) THEN
                           OK1 = .FALSE.
                        END IF
                     END DO
                  END DO
               ELSE
                  DO J = 1, N
                     DO I = 1, J
                        IF ( A(I,J).NE.ASAV(I,J) ) THEN
                           OK1 = .FALSE.
                        END IF
                     END DO
                  END DO
               END IF
*
               NRUN = NRUN + 1
*
               SRNAMT = 'ZTRTTP'
               CALL ZTRTTP( UPLO, N, A, LDA, AP, INFO )
*
               SRNAMT = 'ZTPTTF'
               CALL ZTPTTF( CFORM, UPLO, N, AP, ARF, INFO )
*
               SRNAMT = 'ZTFTTR'
               CALL ZTFTTR( CFORM, UPLO, N, ARF, ASAV, LDA, INFO )
*
               OK2 = .TRUE.
               IF ( LOWER ) THEN
                  DO J = 1, N
                     DO I = J, N
                        IF ( A(I,J).NE.ASAV(I,J) ) THEN
                           OK2 = .FALSE.
                        END IF
                     END DO
                  END DO
               ELSE
                  DO J = 1, N
                     DO I = 1, J
                        IF ( A(I,J).NE.ASAV(I,J) ) THEN
                           OK2 = .FALSE.
                        END IF
                     END DO
                  END DO
               END IF
*
               IF (( .NOT.OK1 ).OR.( .NOT.OK2 )) THEN
                  IF( NERRS.EQ.0 ) THEN
                     WRITE( NOUT, * )
                     WRITE( NOUT, FMT = 9999 )
                  END IF
                  WRITE( NOUT, FMT = 9998 ) N, UPLO, CFORM
                  NERRS = NERRS + 1
               END IF
*
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
*
*     Print a summary of the results.
*
      IF ( NERRS.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9997 ) NRUN
      ELSE
         WRITE( NOUT, FMT = 9996 ) NERRS, NRUN
      END IF
*
 9999 FORMAT( 1X, ' *** Error(s) while testing the RFP conversion',
     +         ' routines ***')
 9998 FORMAT( 1X, '     Error in RFP,conversion routines N=',I5,
     +        ' UPLO=''', A1, ''', FORM =''',A1,'''')
 9997 FORMAT( 1X, 'All tests for the RFP conversion routines passed (',
     +        I5,' tests run)')
 9996 FORMAT( 1X, 'RFP conversion routines:',I5,' out of ',I5,
     +        ' error message recorded')
*
      RETURN
*
*     End of ZDRVRF2
*
      END
