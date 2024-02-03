      int              FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      int                ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN, REAL
*     ..
*     .. External Functions ..
      int                IEEECK, IPARAM2STAGE
      EXTERNAL           IEEECK, IPARAM2STAGE
*     ..
*     .. Arrays in Common ..
      int                IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Save statement ..
      SAVE               / CLAENV /
*     ..
*     .. Executable Statements ..
*
      IF( ISPEC.GE.1 .AND. ISPEC.LE.5 ) THEN
*
*        Return a value from the common block.
*
         ILAENV = IPARMS( ISPEC )
*
      ELSE IF( ISPEC.EQ.6 ) THEN
*
*        Compute SVD crossover point.
*
         ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
*
      ELSE IF( ISPEC.GE.7 .AND. ISPEC.LE.9 ) THEN
*
*        Return a value from the common block.
*
         ILAENV = IPARMS( ISPEC )
*
      ELSE IF( ISPEC.EQ.10 ) THEN
*
*        IEEE NaN arithmetic can be trusted not to trap
*
C        ILAENV = 0
         ILAENV = 1
         IF( ILAENV.EQ.1 ) THEN
            ILAENV = IEEECK( 1, 0.0, 1.0 )
         END IF
*
      ELSE IF( ISPEC.EQ.11 ) THEN
*
*        Infinity arithmetic can be trusted not to trap
*
C        ILAENV = 0
         ILAENV = 1
         IF( ILAENV.EQ.1 ) THEN
            ILAENV = IEEECK( 0, 0.0, 1.0 )
         END IF
*
      ELSE IF(( ISPEC.GE.12 ) .AND. (ISPEC.LE.16)) THEN
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
*
         ILAENV = IPARMS( ISPEC )
*         WRITE(*,*) 'ISPEC = ',ISPEC,' ILAENV =',ILAENV
*         ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
      ELSE IF(( ISPEC.GE.17 ) .AND. (ISPEC.LE.21)) THEN
*
*     17 <= ISPEC <= 21: 2stage eigenvalues SVD routines.
*
         IF( ISPEC.EQ.17 ) THEN
             ILAENV = IPARMS( 1 )
         ELSE
             ILAENV = IPARAM2STAGE( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
         ENDIF
*
      ELSE
*
*        Invalid value for ISPEC
*
         ILAENV = -1
      END IF
*
      RETURN
*
*     End of ILAENV
*
      END
      int     FUNCTION ILAENV2STAGE( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      int                ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Local variables ..
      int                IISPEC
*     .. External Functions ..
      int                IPARAM2STAGE
      EXTERNAL           IPARAM2STAGE
*     ..
*     .. Arrays in Common ..
      int                IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Save statement ..
      SAVE               / CLAENV /
*     ..
*     .. Executable Statements ..
*
      IF(( ISPEC.GE.1 ) .AND. (ISPEC.LE.5)) THEN
*
*     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.
*
         IF( ISPEC.EQ.1 ) THEN
             ILAENV2STAGE = IPARMS( 1 )
         ELSE
             IISPEC = 16 + ISPEC
             ILAENV2STAGE = IPARAM2STAGE( IISPEC, NAME, OPTS, N1, N2, N3, N4 )
         ENDIF
*
      ELSE
*
*        Invalid value for ISPEC
*
         ILAENV2STAGE = -1
      END IF
*
      RETURN
*
*     End of ILAENV2STAGE
*
      END
      int     FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
      int                INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16 )
      int                NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 11, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Scalar Arguments ..
      int                IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*     ..
*     .. Local Scalars ..
      int                NH, NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) NS = 4          IF( NH.GE.60 ) NS = 10          IF( NH.GE.150 ) NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )          IF( NH.GE.590 ) NS = 64          IF( NH.GE.3000 ) NS = 128          IF( NH.GE.6000 ) NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to LAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
         IPARMQ = 0
         IF( NS.GE.KACMIN ) IPARMQ = 1          IF( NS.GE.K22MIN ) IPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
