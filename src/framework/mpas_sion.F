module mpas_sion

   use mpas_dmpar, only : mpas_dmpar_max_int8
   use mpas_derived_types
   use mpas_typecast
   use mpas_log
#ifdef SIONLIB
   use sion_f90
   use sion_f90_mpi
#endif

   implicit none

! Enable debugging information
!define SIONLIB_DEBUG


! Enables prepending variable ID and variable size to the data block for
! reading and writing - highly inefficient, but for debugging etc;
! reads and writes must use consistent setting of this property.
! To guarantee this, an appropriate entry is added to the metadata.
#define SIONLIB_PREPEND_VARIABLE_PROPERTIES

   ! The chunksize for metadata is currently fixed at 40kB, which is more
   ! than twice what is needed for MPAS release 4.0 to store all metadata
   ! for restart files (largest set of metadata). Should in newer versions
   ! variable attributes be added, then this must be checked again.
   ! The chunksize for data files is only used for reading (and will be
   ! overwritten by the SIONlib open routine), for writing the correct
   ! size of the chunksizes is calculated on the fly.
   integer (kind=I8KIND), private :: SIONLIB_metadata_buffer_size = 81920 ! 80kB
   integer (kind=I8KIND), private :: SIONLIB_chunksize_data = 1073741824  ! 1GB, will be overwritten for read and write
   integer,   private :: SIONLIB_fsblksize = -1               ! autodetect

   ! Likewise, the number of sionlib files is controlled by the
   ! namelist parameter config_sionlib_numfiles, which defaults to 1
   integer, private :: SIONLIB_numfiles = 1                   ! can be overwritten by config_sionlib_numfiles

   ! These hard-coded integer values are matching those in PIO 1.9.23
   integer, parameter :: SIONLIB_noerr         = 0            ! PIO_noerr
   integer, parameter :: SIONLIB_err           = -53          ! PIO_RETURN_ERROR
   integer, parameter :: SIONlib_global        = 0            ! PIO_global

   integer, parameter :: HugeStrKIND = 1000

contains

   subroutine MPAS_SION_init(dminfo, io_numfiles, ierr)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: io_numfiles
      integer, intent(out) :: ierr

#ifdef SIONLIB
#ifdef SIONLIB_DEBUG
      call mpas_log_write(' Called MPAS_SION_init')
#endif

      if (io_numfiles < 1 .or. io_numfiles > dminfo % nprocs) then
         call mpas_log_write('Invalid value for config_sionlib_numfiles, must be between 1 and $i', messageType=MPAS_LOG_ERR, intArgs=(/ dminfo % nprocs/))
         SIONLIB_numfiles = 1
         ierr = SIONLIB_err
      else
         SIONLIB_numfiles = io_numfiles
         ierr = SIONLIB_noerr
      end if
#endif

   end subroutine MPAS_SION_init

   subroutine MPAS_SION_open(handle, filename, dminfo, mode, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: filename
      type (dm_info), intent(inout) :: dminfo
      integer, intent(in) :: mode
      integer, intent(out) :: ierr

#ifdef SIONLIB
      ierr = SIONLIB_noerr

      ! Open a SIONlib file. For reading open the file directly. For
      ! writing, defer the opening since the required information on
      ! the storage per task (chunksize) will be calculated on the fly
      ! and only available when MPAS_SION_enddef is called.
      if (handle % sidd /= -1) then
         call mpas_log_write('ERROR in MPAS_SION_open: SIONlib file ' // trim(filename) // ' is already opened', messageType=MPAS_LOG_ERR)
         ierr = MPAS_IO_ERR_OPEN
         return
      end if
      handle % data_mode = .false.
      handle % filename = trim(filename)
      ! Allocate metadata buffer and set cursor to beginning of buffer
      allocate(handle % metadata_buffer(SIONLIB_metadata_buffer_size))
      handle % metadata_cursor = 1
      ! Account for the metadata buffer being written to the beginning
      ! of the file when recording the storage positions of variables
      handle % data_cursor = SIONLIB_metadata_buffer_size

      if (mode == MPAS_IO_READ) then
         call MPAS_SION_open_file(handle % sidd, trim(filename), &
                                 dminfo % comm, dminfo % nullcomm, &
                                 dminfo % my_proc_id, mode, &
                                 SIONLIB_chunksize_data, ierr)
         handle % start_cursor = FSION_GET_POSITION(handle % sidd)
      end if
#endif

   end subroutine MPAS_SION_open

   subroutine MPAS_SION_enddef(handle, filename, dminfo, mode, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: filename
      type (dm_info), intent(inout) :: dminfo
      integer, intent(in) :: mode
      integer, intent(out) :: ierr

      integer (kind=I8KIND) :: position, count
      integer (kind=I8KIND) :: chunksize_max

#ifdef SIONLIB
      if (handle % sidd /= -1) then
         call mpas_log_write('ERROR in MPAS_SION_open: SIONlib file ' // trim(filename) // ' is already opened', messageType=MPAS_LOG_ERR)
         ierr = MPAS_IO_ERR_OPEN
         return
      end if

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      ! This will be tested during reading for consistency checking
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'dbg')
#endif

      ! Gather global maximum of all chunksizes stored in form of current data_cursor
      call mpas_dmpar_max_int8(dminfo, handle % data_cursor, chunksize_max)
      ! Adjust to file system block size
      chunksize_max = ((chunksize_max/SIONLIB_fsblksize)+1)*SIONLIB_fsblksize
#ifdef SIONLIB_DEBUG
      call mpas_log_write('MPAS_SION_enddef -- maximum chunksize is $i', intArgs=(/ int(chunksize_max) /))
#endif
      ! Open the data file for writing using the calculated chunksize
      call MPAS_SION_open_file(handle % sidd, trim(filename), &
                              dminfo % comm, dminfo % nullcomm, &
                              dminfo % my_proc_id, mode, &
                              chunksize_max, ierr)
      if (ierr /= MPAS_IO_NOERR) return
      handle % data_mode = .true.
      handle % start_cursor = FSION_GET_POSITION(handle % sidd)

      ! All metadata collected, write the metadata buffer to the file
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'cur')
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, handle % start_cursor)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'eof')
#ifdef SIONLIB_DEBUG
      call mpas_log_write('MPAS_SION_enddef -- metadata cursor is $i', intArgs=(/ handle % metadata_cursor /) )
#endif
      if (handle % metadata_cursor >= SIONLIB_metadata_buffer_size) then
         call mpas_log_write('ERROR - buffer overflow for SIONlib metadata buffer', messageType=MPAS_LOG_ERR)
         ierr = SIONLIB_err
         return
      end if
      call FSION_COLL_FWRITE_MPI(handle % metadata_buffer(1), int(1,8), SIONLIB_metadata_buffer_size, handle % sidd, count)
      if (count /= SIONLIB_metadata_buffer_size) then
         ierr = SIONLIB_err
         return
      end if
#ifdef SIONLIB_DEBUG
      call mpas_log_write('MPAS_SION_enddef -- wrote $i bytes of SIONlib metadata to disk, new position is now $i', &
                          intArgs=(/ int(count), int(FSION_GET_POSITION(handle % sidd)) /) )
#endif
#endif

   end subroutine MPAS_SION_enddef

   subroutine MPAS_SION_open_file(sid, filename, comm, nullcomm, proc_id, mode, chunksize, ierr)

      implicit none

      integer, intent(out) :: sid
      character (len=*), intent(in) :: filename
      integer, intent(in) :: comm
      integer, intent(inout) :: nullcomm
      integer, intent(in) :: proc_id
      integer (kind=I8KIND), intent(inout) :: chunksize
      integer, intent(in) :: mode

      integer, intent(out) :: ierr

#ifdef SIONLIB
      ! Maximum length of filename as for filename in mpas_io.F
      character(len=1024) :: newfname
      integer :: f_endian, s_endian
      character(len=2) :: cmode
      integer :: sionlib_numfiles_local

#ifdef SIONLIB_DEBUG
      call mpas_log_write(' Called MPAS_SION_open_file')
#endif

      ! Required so that an open-for-reading is not overwriting the user-defined value
      sionlib_numfiles_local = SIONLIB_numfiles
      if (mode == MPAS_IO_WRITE) then
          cmode = 'wb'
      else
          cmode = 'rb'
      end if

      ! Parallel open file, if successful sid >= 0
      call FSION_PAROPEN_MPI(trim(filename), cmode, sionlib_numfiles_local, &
                             comm, nullcomm, chunksize, SIONLIB_fsblksize, &
                             proc_id, newfname, sid)
      if (sid .eq. -1) then
          ierr = MPAS_IO_ERR_OPEN
          return
      else
          ierr = MPAS_IO_NOERR
      end if

      ! Check that file endianness is identical to system endianness
      call FSION_GET_FILE_ENDIANNESS(sid, f_endian)
      call FSION_GET_ENDIANESS(s_endian)
      if (f_endian .ne. s_endian) then
         ierr = MPAS_IO_ERR_ENDIANNESS
         return
      end if

#ifdef SIONLIB_DEBUG
      call mpas_log_write('Opening file ' // trim(filename) // ' in mode "' // trim(cmode) // '" as set of $i files', intArgs=(/ SIONLIB_numfiles /))
#endif
#endif

   end subroutine MPAS_SION_open_file

   subroutine MPAS_SION_def_dim(handle, dimname, dimid, dimsize, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: dimname
      integer, intent(in) :: dimid
      integer, intent(in) :: dimsize
      integer, intent(out) :: ierr

#ifdef SIONLIB
      integer :: cursor_start, checksum

      cursor_start = handle % metadata_cursor
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'dim')
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, dimid)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, len(trim(dimname)))
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, trim(dimname))
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, dimsize)
      checksum = 3 + kind(dimid) + kind(len(trim(dimname))) + len(trim(dimname)) + kind(dimsize)

#ifdef SIONLIB_DEBUG
      call mpas_log_write('After encoding dimension ' // trim(dimname) // ', cursor = $i', intArgs=(/ handle % metadata_cursor /) )
#endif
      if (handle % metadata_cursor - cursor_start == checksum) then
          ierr = SIONLIB_noerr
      else
          ierr = SIONLIB_err
      end if
#endif

   end subroutine MPAS_SION_def_dim

   subroutine MPAS_SION_def_var(handle, varname, varid, ndims, dimids, vartype, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: varname
      integer, intent(in) :: varid
      integer, intent(in) :: ndims
      integer, dimension(1:ndims), intent(in) :: dimids
      integer, intent(in) :: vartype
      integer, intent(out) :: ierr

#ifdef SIONLIB
      integer :: idim, cursor_start, checksum

      cursor_start = handle % metadata_cursor
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'var')
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, varid)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, len(trim(varname)))
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, trim(varname))
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, ndims)
      do idim=1,ndims
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, dimids(idim))
      end do
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, vartype)
      checksum = 3 + kind(varid) + kind(len(trim(varname))) + len(trim(varname)) + kind(ndims) &
               & + ndims*kind(ndims) + kind(vartype)

#ifdef SIONLIB_DEBUG
      call mpas_log_write('After encoding variable ' // trim(varname) // ', cursor = $i', intArgs=(/ handle % metadata_cursor /) )
#endif
      if (handle % metadata_cursor - cursor_start == checksum) then
          ierr = SIONLIB_noerr
      else
          ierr = SIONLIB_err
      end if
#endif

   end subroutine MPAS_SION_def_var

   subroutine MPAS_SION_initdecomp(decomp, dimlist, ierr)

      implicit none

      type (decomplist_type), pointer, intent(inout) :: decomp
      integer, dimension(:), pointer, intent(in) :: dimlist
      integer, intent(out) :: ierr

#ifdef SIONLIB
      ierr = SIONLIB_noerr

      ! Now calculate the chunksize for this given decomposition type, using
      ! the information on the field type and the dimensions of the decomp handle
      decomp % decomphandle % sionlib_fieldsize = product(dimlist)
      call MPAS_SION_chunksize(decomp % decomphandle % field_type, &
                               decomp % decomphandle % sionlib_fieldsize, &
                               decomp % decomphandle % sionlib_chunksize, ierr)
      if (ierr /= SIONLIB_noerr) return

      decomp % decomphandle % sionlib_type = .true.
#endif

      !DIR$ ATTRIBUTES INLINE :: MPAS_SION_initdecomp
   end subroutine MPAS_SION_initdecomp

   subroutine MPAS_SION_chunksize(vartype, varsize, chunksize, ierr)

      implicit none

      integer, intent(in) :: vartype
      integer (kind=I8KIND), intent(in) :: varsize
      integer (kind=I8KIND), intent(out) :: chunksize
      integer, intent(out) :: ierr

#ifdef SIONLIB
      ierr = SIONLIB_noerr

      if (vartype == MPAS_IO_REAL) then
         chunksize = RKIND * varsize
      else if (vartype == MPAS_IO_DOUBLE) then
         chunksize = R8KIND * varsize
      else if (vartype == MPAS_IO_INT) then
         chunksize = kind(SIONLIB_noerr) * varsize
      else if (vartype == MPAS_IO_CHAR) then
         chunksize = varsize
      else
         call mpas_log_write('ERROR in MPAS_SION_def_chunksize, fieldtype $i not supported', messageType=MPAS_LOG_ERR, intArgs=(/ vartype /) )
         chunksize = 0
         ierr = SIONLIB_err
      end if

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      chunksize = chunksize + &
                  kind(vartype)*2 + & ! varid and vartype
                  kind(varsize)*1     ! varsize
#endif
#endif

   end subroutine MPAS_SION_chunksize

   subroutine MPAS_SION_store_position_non_decomp(handle, field_cursor, fieldsize, chunksize, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      type (fieldlist_type), pointer, intent(inout) :: field_cursor
      integer (kind=I8KIND), intent(in) :: fieldsize
      integer (kind=I8KIND), intent(in) :: chunksize
      integer, intent(inout) :: ierr

#ifdef SIONLIB
      integer :: cursor_start
      integer :: checksum

#ifdef SIONLIB_DEBUG
      call mpas_log_write('Storing position $i for non-decomposed variable $i / ' // &
         trim(field_cursor % fieldhandle % fieldname) // ' of size $i, advance by $i', &
         intArgs=(/ int(handle % data_cursor), field_cursor % fieldhandle % fieldid, int(fieldsize), int(chunksize) /) )
#endif
      field_cursor % fieldhandle % sionlib_position  = handle % data_cursor
      field_cursor % fieldhandle % sionlib_fieldsize = fieldsize
      field_cursor % fieldhandle % sionlib_chunksize = chunksize
      handle % data_cursor = handle % data_cursor + chunksize
      ! Encode metadata
      cursor_start = handle % metadata_cursor
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'pos')
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, field_cursor % fieldhandle % fieldid)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, field_cursor % fieldhandle % sionlib_position)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, fieldsize)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, chunksize)
      checksum = 3 + kind(field_cursor % fieldhandle % fieldid) + &
               & kind(field_cursor % fieldhandle % sionlib_position) + &
               & kind(fieldsize) + &
               & kind(chunksize)
#ifdef SIONLIB_DEBUG
      call mpas_log_write('After encoding position for non-decomposed varid $i, cursor = $i', &
                     intArgs=(/ field_cursor % fieldhandle % fieldid, handle % metadata_cursor /) )
#endif
      if (handle % metadata_cursor - cursor_start == checksum) then
          ierr = SIONLIB_noerr
      else
          ierr = SIONLIB_err
          return
      end if
#endif
   end subroutine MPAS_SION_store_position_non_decomp

   subroutine MPAS_SION_store_position_decomp(handle, field_cursor, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      type (fieldlist_type), pointer, intent(inout) :: field_cursor
      integer, intent(inout) :: ierr

#ifdef SIONLIB
      integer :: cursor_start
      integer :: checksum

#ifdef SIONLIB_DEBUG
      call mpas_log_write('Storing position $i for decomposed variable $i / ' // &
         trim(field_cursor % fieldhandle % fieldname) // ' of size $i, advance by $i', &
         intArgs=(/ int(handle % data_cursor), field_cursor % fieldhandle % fieldid, &
                    int(field_cursor % fieldhandle % decomp % sionlib_fieldsize), &
                    int(field_cursor % fieldhandle % decomp % sionlib_chunksize) /) )
#endif
      field_cursor % fieldhandle % sionlib_position  = handle % data_cursor
      field_cursor % fieldhandle % sionlib_fieldsize = field_cursor % fieldhandle % decomp % sionlib_fieldsize
      field_cursor % fieldhandle % sionlib_chunksize = field_cursor % fieldhandle % decomp % sionlib_chunksize
      handle % data_cursor = handle % data_cursor + field_cursor % fieldhandle % decomp % sionlib_chunksize
      ! Encode metadata
      cursor_start = handle % metadata_cursor
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'pos')
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, field_cursor % fieldhandle % fieldid)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, field_cursor % fieldhandle % sionlib_position)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, field_cursor % fieldhandle % decomp % sionlib_fieldsize)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, field_cursor % fieldhandle % decomp % sionlib_chunksize)
      checksum = 3 + kind(field_cursor % fieldhandle % fieldid) + &
               & kind(field_cursor % fieldhandle % sionlib_position)  + &
               & kind(field_cursor % fieldhandle % decomp % sionlib_fieldsize) + &
               & kind(field_cursor % fieldhandle % decomp % sionlib_chunksize)
#ifdef SIONLIB_DEBUG
      call mpas_log_write('After encoding position for decomposed varid $i, cursor = $i', &
                     intArgs=(/ field_cursor % fieldhandle % fieldid, handle % metadata_cursor /) )
#endif
      if (handle % metadata_cursor - cursor_start == checksum) then
          ierr = SIONLIB_noerr
      else
          ierr = SIONLIB_err
          return
      end if
#endif
   end subroutine MPAS_SION_store_position_decomp

   subroutine MPAS_SION_write_array(handle, field_cursor, &
      & realval, realarray1d, realarray2d, realarray3d, realarray4d, realarray5d, &
      & intval, intarray1d, intarray2d, intarray3d, intarray4d, intarray5d, &
      & singleval, singlearray1d, singlearray2d, singlearray3d, singlearray4d, singlearray5d, &
      & doubleval, doublearray1d, doublearray2d, doublearray3d, doublearray4d, doublearray5d, &
      & charval, chararray1d, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      !integer, intent(in) :: varid
      type (fieldlist_type), pointer, intent(inout) :: field_cursor
      real (kind=RKIND),  intent(in), optional :: realval
      real (kind=RKIND),  dimension(:), intent(in), optional :: realarray1d
      real (kind=RKIND),  dimension(:,:), intent(in), optional :: realarray2d
      real (kind=RKIND),  dimension(:,:,:), intent(in), optional :: realarray3d
      real (kind=RKIND),  dimension(:,:,:,:), intent(in), optional :: realarray4d
      real (kind=RKIND),  dimension(:,:,:,:,:), intent(in), optional :: realarray5d
      integer,            intent(in), optional :: intval
      integer,            dimension(:), intent(in), optional :: intarray1d
      integer,            dimension(:,:), intent(in), optional :: intarray2d
      integer,            dimension(:,:,:), intent(in), optional :: intarray3d
      integer,            dimension(:,:,:,:), intent(in), optional :: intarray4d
      integer,            dimension(:,:,:,:,:), intent(in), optional :: intarray5d
      real (kind=R4KIND), intent(in), optional :: singleval
      real (kind=R4KIND), dimension(:), intent(in), optional :: singlearray1d
      real (kind=R4KIND), dimension(:,:), intent(in), optional :: singlearray2d
      real (kind=R4KIND), dimension(:,:,:), intent(in), optional :: singlearray3d
      real (kind=R4KIND), dimension(:,:,:,:), intent(in), optional :: singlearray4d
      real (kind=R4KIND), dimension(:,:,:,:,:), intent(in), optional :: singlearray5d
      real (kind=R8KIND), intent(in), optional :: doubleval
      real (kind=R8KIND), dimension(:), intent(in), optional :: doublearray1d
      real (kind=R8KIND), dimension(:,:), intent(in), optional :: doublearray2d
      real (kind=R8KIND), dimension(:,:,:), intent(in), optional :: doublearray3d
      real (kind=R8KIND), dimension(:,:,:,:), intent(in), optional :: doublearray4d
      real (kind=R8KIND), dimension(:,:,:,:,:), intent(in), optional :: doublearray5d
      character (len=*),  intent(in), optional :: charval
      character (len=*),  dimension(:), intent(in), optional :: chararray1d

      integer, intent(out) :: ierr

#ifdef SIONLIB
      integer :: checksum
      integer (kind=I8KIND) :: count, position
      character (len=StrKIND) :: padded_charval
      character (len=HugeStrKIND) :: padded_charval_huge
      integer :: ichar

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      integer :: varid
      integer :: vartype
      integer (kind=I8KIND) :: varsize
#endif

      ! Get current position in file and check that it matches
      ! the calculated position stored in the metadata
      position = FSION_GET_POSITION(handle % sidd)

      if (field_cursor % fieldhandle % sionlib_position /= position - handle % start_cursor) then
         call mpas_log_write('ERROR: calculated position for writing to file differs from actual position: $i vs $i', &
                        messageType=MPAS_LOG_ERR, intArgs=(/ int(field_cursor % fieldhandle % sionlib_position), int(position - handle % start_cursor) /) )
         ierr = SIONLIB_err
         return
      end if

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      varid   = field_cursor % fieldhandle % fieldid
      vartype = field_cursor % fieldhandle % field_type
      varsize = field_cursor % fieldhandle % sionlib_fieldsize
      call FSION_COLL_FWRITE_MPI(varid,   int(kind(varid),8),   int(1,8), handle % sidd, count)
      call FSION_COLL_FWRITE_MPI(vartype, int(kind(vartype),8), int(1,8), handle % sidd, count)
      call FSION_COLL_FWRITE_MPI(varsize, int(kind(varsize),8), int(1,8), handle % sidd, count)
#endif

      if (present(realval)) then
         checksum = 1
         call FSION_COLL_FWRITE_MPI(realval, int(kind(realval),8), int(1,8), handle % sidd, count)
      else if (present(realarray1d)) then
         checksum = size(realarray1d)
         call FSION_COLL_FWRITE_MPI(realarray1d(1), int(kind(realarray1d(1)),8), int(size(realarray1d),8), handle % sidd, count)
      else if (present(realarray2d)) then
         checksum = size(realarray2d)
         call FSION_COLL_FWRITE_MPI(realarray2d(1,1), int(kind(realarray2d(1,1)),8), int(size(realarray2d),8), handle % sidd, count)
      else if (present(realarray3d)) then
         checksum = size(realarray3d)
         call FSION_COLL_FWRITE_MPI(realarray3d(1,1,1), int(kind(realarray3d(1,1,1)),8), int(size(realarray3d),8), handle % sidd, count)
      else if (present(realarray4d)) then
         checksum = size(realarray4d)
         call FSION_COLL_FWRITE_MPI(realarray4d(1,1,1,1), int(kind(realarray4d(1,1,1,1)),8), int(size(realarray4d),8), handle % sidd, count)
      else if (present(realarray5d)) then
         checksum = size(realarray5d)
         call FSION_COLL_FWRITE_MPI(realarray5d(1,1,1,1,1), int(kind(realarray5d(1,1,1,1,1)),8), int(size(realarray5d),8), handle % sidd, count)
      else if (present(intval)) then
         checksum = 1
         call FSION_COLL_FWRITE_MPI(intval, int(kind(intval),8), int(1,8), handle % sidd, count)
      else if (present(intarray1d)) then
         checksum = size(intarray1d)
         call FSION_COLL_FWRITE_MPI(intarray1d(1), int(kind(intarray1d(1)),8), int(size(intarray1d),8), handle % sidd, count)
      else if (present(intarray2d)) then
         checksum = size(intarray2d)
         call FSION_COLL_FWRITE_MPI(intarray2d(1,1), int(kind(intarray2d(1,1)),8), int(size(intarray2d),8), handle % sidd, count)
      else if (present(intarray3d)) then
         checksum = size(intarray3d)
         call FSION_COLL_FWRITE_MPI(intarray3d(1,1,1), int(kind(intarray3d(1,1,1)),8), int(size(intarray3d),8), handle % sidd, count)
      else if (present(intarray4d)) then
         checksum = size(intarray4d)
         call FSION_COLL_FWRITE_MPI(intarray4d(1,1,1,1), int(kind(intarray4d(1,1,1,1)),8), int(size(intarray4d),8), handle % sidd, count)
      else if (present(intarray5d)) then
         checksum = size(intarray5d)
         call FSION_COLL_FWRITE_MPI(intarray5d(1,1,1,1,1), int(kind(intarray5d(1,1,1,1,1)),8), int(size(intarray5d),8), handle % sidd, count)
      else if (present(singleval)) then
         checksum = 1
         call FSION_COLL_FWRITE_MPI(singleval, int(kind(singleval),8), int(1,8), handle % sidd, count)
      else if (present(singlearray1d)) then
         checksum = size(singlearray1d)
         call FSION_COLL_FWRITE_MPI(singlearray1d(1), int(kind(singlearray1d(1)),8), int(size(singlearray1d),8), handle % sidd, count)
      else if (present(singlearray2d)) then
         checksum = size(singlearray2d)
         call FSION_COLL_FWRITE_MPI(singlearray2d(1,1), int(kind(singlearray2d(1,1)),8), int(size(singlearray2d),8), handle % sidd, count)
      else if (present(singlearray3d)) then
         checksum = size(singlearray3d)
         call FSION_COLL_FWRITE_MPI(singlearray3d(1,1,1), int(kind(singlearray3d(1,1,1)),8), int(size(singlearray3d),8), handle % sidd, count)
      else if (present(singlearray4d)) then
         checksum = size(singlearray4d)
         call FSION_COLL_FWRITE_MPI(singlearray4d(1,1,1,1), int(kind(singlearray4d(1,1,1,1)),8), int(size(singlearray4d),8), handle % sidd, count)
      else if (present(singlearray5d)) then
         checksum = size(singlearray5d)
         call FSION_COLL_FWRITE_MPI(singlearray5d(1,1,1,1,1), int(kind(singlearray5d(1,1,1,1,1)),8), int(size(singlearray5d),8), handle % sidd, count)
      else if (present(doubleval)) then
         checksum = 1
         call FSION_COLL_FWRITE_MPI(doubleval, int(kind(doubleval),8), int(1,8), handle % sidd, count)
      else if (present(doublearray1d)) then
         checksum = size(doublearray1d)
         call FSION_COLL_FWRITE_MPI(doublearray1d(1), int(kind(doublearray1d(1)),8), int(size(doublearray1d),8), handle % sidd, count)
      else if (present(doublearray2d)) then
         checksum = size(doublearray2d)
         call FSION_COLL_FWRITE_MPI(doublearray2d(1,1), int(kind(doublearray2d(1,1)),8), int(size(doublearray2d),8), handle % sidd, count)
      else if (present(doublearray3d)) then
         checksum = size(doublearray3d)
         call FSION_COLL_FWRITE_MPI(doublearray3d(1,1,1), int(kind(doublearray3d(1,1,1)),8), int(size(doublearray3d),8), handle % sidd, count)
      else if (present(doublearray4d)) then
         checksum = size(doublearray4d)
         call FSION_COLL_FWRITE_MPI(doublearray4d(1,1,1,1), int(kind(doublearray4d(1,1,1,1)),8), int(size(doublearray4d),8), handle % sidd, count)
      else if (present(doublearray5d)) then
         checksum = size(doublearray5d)
         call FSION_COLL_FWRITE_MPI(doublearray5d(1,1,1,1,1), int(kind(doublearray5d(1,1,1,1,1)),8), int(size(doublearray5d),8), handle % sidd, count)
      else if (present(charval)) then
         padded_charval(1:len(trim(charval))) = trim(charval)
         padded_charval(len(trim(charval))+1:) = ' '
         checksum = len(padded_charval)
         call FSION_COLL_FWRITE_MPI(padded_charval, int(1,8), int(len(padded_charval),8), handle % sidd, count)
      else if (present(chararray1d)) then
         checksum = size(chararray1d)*len(padded_charval)
         if (checksum > len(padded_charval_huge)) then
             call mpas_log_write('MPAS_SION_write_darray: error, buffer overflow for writing chararray1d from file for variable $i - increase HugeStrKIND!', &
                            messageType=MPAS_LOG_ERR, intArgs=(/ field_cursor % fieldhandle % fieldid /) )
             ierr = SIONLIB_err
             return
         end if
         ! Encode character array as one long string in padded_charval_huge and write to file
         do ichar=1,size(chararray1d)
            padded_charval(1:len(trim(charval))) = trim(chararray1d(ichar))
            padded_charval(len(trim(chararray1d(ichar)))+1:) = ' '
            padded_charval_huge((ichar-1)*len(padded_charval)+1:ichar*len(padded_charval)) = padded_charval
         end do
         call FSION_COLL_FWRITE_MPI(padded_charval_huge, int(1,8), int(len(padded_charval_huge),8), handle % sidd, count)
         ! DH*
         call mpas_log_write('Attention, writing charArray1d variables has not been tested yet!', messageType=MPAS_LOG_WARN)
         ! *DH
      else
          call mpas_log_write('MPAS_SION_write_darray: no values specified for variable $i', messageType=MPAS_LOG_ERR, intArgs=(/ field_cursor % fieldhandle % fieldid /) )
          ierr = SIONLIB_err
          return
      end if

      if (count == checksum) then
          ierr = SIONLIB_noerr
      else
          ierr = SIONLIB_err
      end if
#endif

   end subroutine MPAS_SION_write_array

   subroutine MPAS_SION_put_att(handle, varid, attname, atttype, charval, intval, intarray1d, &
                                realval, realarray1d, singleval, singlearray1d, doubleval, &
                                doublearray1d, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(in) :: varid
      character (len=*), intent(in) :: attname
      integer, intent(in) :: atttype
      character (len=*), intent(in), optional :: charval
      integer, intent(in), optional :: intval
      integer, dimension(:), intent(in), optional :: intarray1d
      real (kind=RKIND), intent(in), optional :: realval
      real (kind=RKIND), dimension(:), intent(in), optional :: realarray1d
      real (kind=R4KIND), intent(in), optional :: singleval
      real (kind=R4KIND), dimension(:), intent(in), optional :: singlearray1d
      real (kind=R8KIND), intent(in), optional :: doubleval
      real (kind=R8KIND), dimension(:), intent(in), optional :: doublearray1d

      integer, intent(out) :: ierr

#ifdef SIONLIB
      integer :: iatt, cursor_start, checksum

      cursor_start = handle % metadata_cursor
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 'att')
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, varid)
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, len(trim(attname)))
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, trim(attname))
      call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, atttype)
      ! Encode native precision for all variables other than single/double
      if (present(singleval) .or. present(singlearray1d)) then
          call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, MPAS_IO_SINGLE_PRECISION)
      else if (present(doubleval) .or. present(doublearray1d)) then
          call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, MPAS_IO_DOUBLE_PRECISION)
      else
          call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, MPAS_IO_NATIVE_PRECISION)
      end if
      checksum = 3 + kind(varid) + kind(len(trim(attname))) + len(trim(attname)) + kind(atttype) + kind(MPAS_IO_NATIVE_PRECISION)

      if (present(intval)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 1)
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, intval)
         checksum = checksum + kind(1) + kind(intval)
      else if (present(intarray1d)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, size(intarray1d))
         do iatt=1,size(intarray1d)
            call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, intarray1d(iatt))
         end do
         checksum = checksum + kind(1) + size(intarray1d)*kind(intarray1d(1))
      else if (present(realval)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 1)
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, realval)
         checksum = checksum + kind(1) + kind(realval)
      else if (present(realarray1d)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, size(realarray1d))
         do iatt=1,size(realarray1d)
            call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, realarray1d(iatt))
         end do
         checksum = checksum + kind(1) + size(realarray1d)*kind(realarray1d(1))
      else if (present(singleval)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 1)
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, singleval)
         checksum = checksum + kind(1) + kind(singleval)
      else if (present(singlearray1d)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, size(singlearray1d))
         do iatt=1,size(singlearray1d)
            call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, singlearray1d(iatt))
         end do
         checksum = checksum + kind(1) + size(singlearray1d)*kind(singlearray1d(1))
     else if (present(doubleval)) then
        call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, 1)
        call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, doubleval)
        checksum = checksum + kind(1) + kind(doubleval)
     else if (present(doublearray1d)) then
        call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, size(doublearray1d))
        do iatt=1,size(doublearray1d)
           call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, doublearray1d(iatt))
        end do
         checksum = checksum + kind(1) + size(doublearray1d)*kind(doublearray1d(1))
      else if (present(charval)) then
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, len(trim(charval)))
         call mpas_encode(handle % metadata_buffer, handle % metadata_cursor, trim(charval))
         checksum = checksum + kind(len(trim(charval))) + len(trim(charval))
      else
          call mpas_log_write('MPAS_SION_put_att_generic: no value specified for attribute ' // trim(attname), messageType=MPAS_LOG_ERR)
          ierr = SIONLIB_err
          return
      end if

#ifdef SIONLIB_DEBUG
      call mpas_log_write('After encoding attribute ' // trim(attname) // ', cursor = $i', intArgs=(/ handle % metadata_cursor /) )
#endif
      if (handle % metadata_cursor - cursor_start == checksum) then
          ierr = SIONLIB_noerr
      else
          ierr = SIONLIB_err
      end if
#endif

   end subroutine MPAS_SION_put_att

   subroutine MPAS_SION_sync(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

#ifdef SIONLIB
      integer :: local_ierr

      ! Sync files to disk
      if (handle % sidd /= -1) then
         call FSION_FLUSH(handle % sidd, local_ierr)
         if ((local_ierr.ne.1) .and. present(ierr)) ierr = MPAS_IO_ERR_FLUSH
      end if
#endif

   end subroutine MPAS_SION_sync

   subroutine MPAS_SION_close(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

#ifdef SIONLIB
      integer :: local_ierr

      if (handle % sidd == -1) then
          ! File already closed, for example during call to sync
#ifdef SIONLIB_DEBUG
          ! DH* this messsage is written after the log files are close? ends up in fort.1
          ! see similar messages in pull request for interval_in, interval_out support
          call mpas_log_write('Ignoring repeated call to close file ' // trim(handle % filename))
#endif
          return
      end if
#ifdef SIONLIB_DEBUG
      call mpas_log_write('Closing file ' // trim(handle % filename))
#endif

      deallocate(handle % metadata_buffer)
      handle % metadata_cursor = 1

      ! FSION_PARCLOSE_MPI returns 1 if close operation was successful
      call FSION_PARCLOSE_MPI(handle % sidd, local_ierr)
      if ((local_ierr.ne.1) .and. present(ierr)) ierr = MPAS_IO_ERR_CLOSE

      ! Reset file handle id to default/uninitialised value
      handle % sidd = -1
#endif

   end subroutine MPAS_SION_close

   subroutine MPAS_SION_read_metadata(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out) :: ierr

#ifdef SIONLIB
      ! Category of metadata, read buffers and general variables
      integer :: sion_eof
      character (len=3) :: category
      integer (kind=I8KIND) :: count
      integer :: clength
      logical :: found
      integer (kind=I8KIND) :: test_start_cursor

      ! Dimensions
      type (dimlist_type), pointer :: dim_cursor
      integer :: dimid
      character (len=1024) :: dimname
      integer :: dimsize
      integer :: idim

      ! Variables
      type (fieldlist_type), pointer :: field_cursor
      integer :: varid
      character (len=1024) :: varname
      integer :: ndims
      integer, dimension(:), allocatable :: dimids
      integer :: vartype
      integer (kind=I8KIND) :: position
      integer (kind=I8KIND) :: varsize
      integer (kind=I8KIND) :: chunksize

      ! Attributes
      type (attlist_type), pointer :: att_cursor
      character (len=1024) :: attname
      integer :: atttype
      integer :: precision
      integer :: attlength
      integer :: intval
      integer, dimension(:), allocatable :: intarray
      real (KIND=RKIND) :: realval
      real (KIND=RKIND), dimension(:), allocatable :: realarray
      real (KIND=R4KIND) :: singleval
      real (KIND=R4KIND), dimension(:), allocatable :: singlearray
      real (KIND=R8KIND) :: doubleval
      real (KIND=R8KIND), dimension(:), allocatable :: doublearray
      character (len=1024) :: charval
      integer :: iatt
      character (len=HugeStrKIND) :: message
#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      logical :: prepend_check
#endif

      ierr = SIONLIB_noerr

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      prepend_check = .false.
#endif

      ! Read entire metadata contents into buffer
      if (handle % metadata_cursor /= 1) then
         call mpas_log_write('ERROR while reading SIONlib metadata - cursor in wrong position', messageType=MPAS_LOG_ERR)
         ierr = SIONLIB_err
         return
      end if

      ! Make sure that the file cursor is at the start of the block
      call MPAS_SION_set_cursor(handle, handle % start_cursor, ierr)
      if (ierr /= SIONLIB_noerr) return

      call FSION_COLL_FREAD_MPI(handle % metadata_buffer(1), int(1,8), SIONLIB_metadata_buffer_size, handle % sidd, count)
      if (count /= SIONLIB_metadata_buffer_size) then
         call mpas_log_write('Error while reading metadata into buffer, countsize does not match: $i vs $i', &
                        messageType=MPAS_LOG_ERR, intArgs=(/ int(count), int(SIONLIB_metadata_buffer_size) /) )
         ierr = SIONLIB_err
         return
      end if

      do
         call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, category)
         if (trim(category) == 'dim') then
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, dimid)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, clength)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, dimname(1:clength))
            dimname(clength+1:) = ' '
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, dimsize)
#ifdef SIONLIB_DEBUG
            call mpas_log_write('Reading dimension data ...')
            call mpas_log_write('clength               = $i', intArgs=(/ clength /) )
            call mpas_log_write('dimid                 = $i', intArgs=(/ dimid /) )
            call mpas_log_write('dimname               = ' // trim(dimname))
            call mpas_log_write('dimsize               = $i', intArgs=(/ dimsize /) )
            call mpas_log_write('is_unlimited_dim      = $l', logicArgs=(/ dimsize == MPAS_IO_UNLIMITED_DIM /) )
#endif
            call MPAS_SION_read_metadata_def_dim(handle, trim(dimname), dimid, dimsize)

         else if (trim(category) == 'att') then
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, varid)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, clength)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, attname(1:clength))
            attname(clength+1:) = ' '
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, atttype)
            ! Decode precision in file, ignored for non-real attributes
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, precision)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, attlength)
            if (atttype == MPAS_ATT_INT) then
               call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, intval)
               call MPAS_SION_read_metadata_def_att(handle, trim(attname), varid, atttype, intval=intval, ierr=ierr)
            else if (atttype == MPAS_ATT_INTA) then
               allocate(intarray(attlength))
               do iatt=1,attlength
                  call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, intarray(iatt))
               end do
               call MPAS_SION_read_metadata_def_att(handle, trim(attname), varid, atttype, intarray=intarray, ierr=ierr)
               deallocate(intarray)
            else if (atttype == MPAS_ATT_REAL) then
               !call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, realval)
               if (precision == MPAS_IO_SINGLE_PRECISION) then
                   call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, singleval)
                   realval = real(singleval, RKIND)
               else if (precision == MPAS_IO_DOUBLE_PRECISION) then
                   call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, doubleval)
                   realval = real(doubleval, RKIND)
               else
                   call mpas_log_write('Invalid attribute precision encountered in file for attribute ' // trim(attname), messageType=MPAS_LOG_ERR)
                   ierr = SIONLIB_err
                   return
               end if
               call MPAS_SION_read_metadata_def_att(handle, trim(attname), varid, atttype, realval=realval, ierr=ierr)
            else if (atttype == MPAS_ATT_INTA) then
               allocate(realarray(attlength))
               if (precision == MPAS_IO_SINGLE_PRECISION) then
                   allocate(singlearray(attlength))
                   do iatt=1,attlength
                      call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, singlearray(iatt))
                   end do
                   realarray(:) = real(singlearray(:), RKIND)
                   deallocate(singlearray)
               else if (precision == MPAS_IO_DOUBLE_PRECISION) then
                  allocate(doublearray(attlength))
                  do iatt=1,attlength
                     call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, doublearray(iatt))
                  end do
                  realarray(:) = real(doublearray(:), RKIND)
                  deallocate(doublearray)
               else
                   call mpas_log_write('Invalid attribute precision encountered in file for attribute ' // trim(attname), messageType=MPAS_LOG_ERR)
                   ierr = SIONLIB_err
                   return
               end if
               call MPAS_SION_read_metadata_def_att(handle, trim(attname), varid, atttype, realarray=realarray, ierr=ierr)
               deallocate(realarray)
            else if (atttype == MPAS_ATT_TEXT) then
               call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, charval(1:attlength))
               charval(attlength+1:) = ' '
               call MPAS_SION_read_metadata_def_att(handle, trim(attname), varid, atttype, charval=charval, ierr=ierr)
            else
               call mpas_log_write('Invalid attribute type encountered in file for attribute ' // trim(attname), messageType=MPAS_LOG_ERR)
               ierr = SIONLIB_err
            end if
#ifdef SIONLIB_DEBUG
            call mpas_log_write('Reading attribute data ...')
            call mpas_log_write('varid          = $i', intArgs=(/ varid /) )
            call mpas_log_write('clength        = $i', intArgs=(/ clength /) )
            call mpas_log_write('attname        = ' // trim(attname))
            call mpas_log_write('atttype        = $i', intArgs=(/ atttype /) )
            call mpas_log_write('precision      = $i', intArgs=(/ precision /) )
            call mpas_log_write('attlength      = $i', intArgs=(/ attlength /) )
            if (atttype == MPAS_ATT_INT) then
               call mpas_log_write('attval         = $i', intArgs=(/ intval /) )
            else if (atttype == MPAS_ATT_INTA) then
               ! DH* intarray args?
               message = 'attval         = '
               do iatt=1,attlength
                  message = message // '$i '
               end do
               message = trim(message)
               call mpas_log_write(message, intArgs=intarray)
               ! *DH
            else if (atttype == MPAS_ATT_REAL) then
               call mpas_log_write('attval         = $r', realArgs=(/ realval /) )
            else if (atttype == MPAS_ATT_INTA) then
               ! DH* realarrayargs?
               message = 'attval         = '
               do iatt=1,attlength
                  message = message // '$r '
               end do
               message = trim(message)
               call mpas_log_write(message, realArgs=realarray)
               ! *DH
            else if (atttype == MPAS_ATT_TEXT) then
               call mpas_log_write('attval         = ' // trim(charval))
            end if
#endif
            if (ierr /= SIONLIB_noerr) exit

         else if (trim(category) == 'var') then
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, varid)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, clength)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, varname(1:clength))
            varname(clength+1:) = ' '
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, ndims)
            if (ndims .gt. 0) then
               allocate(dimids(1:ndims))
               do idim=1,ndims
                  call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, dimids(idim))
               end do
            end if
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, vartype)
#ifdef SIONLIB_DEBUG
            call mpas_log_write('Reading variable metadata ... ')
            call mpas_log_write('varid                    = $i', intArgs=(/ varid /) )
            call mpas_log_write('clength                  = $i', intArgs=(/ clength /) )
            call mpas_log_write('varname                  = ' // trim(varname))
            call mpas_log_write('ndims                    = $i', intArgs=(/ ndims /) )
            ! DH* intarrayargs?
            if (allocated(dimids)) then
               message = 'attval         = '
               do iatt=1,size(dimids)
                  message = message // '$i '
               end do
               message = trim(message)
               call mpas_log_write(message, intArgs=dimids)
            end if
            ! *DH
            call mpas_log_write('vartype                  = $i', intArgs=(/ vartype /) )
#endif
            call MPAS_SION_read_metadata_def_var(handle, varname, varid, vartype, ndims, dimids, ierr)
            if (allocated(dimids)) deallocate(dimids)
            if (ierr /= SIONLIB_noerr) exit

         else if (trim(category) == 'pos') then
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, varid)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, position)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, varsize)
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, chunksize)

#ifdef SIONLIB_DEBUG
            call mpas_log_write('Reading variable position metadata ... ')
            call mpas_log_write('varid                    = $i', intArgs=(/ varid /) )
            call mpas_log_write('position                 = $i', intArgs=(/ int(position) /) )
            call mpas_log_write('varsize                  = $i', intArgs=(/ int(varsize) /) )
            call mpas_log_write('chunksize                = $i', intArgs=(/ int(chunksize) /) )
#endif
            call MPAS_SION_read_metadata_set_pos(handle, varid, position, varsize, chunksize, ierr)
            if (ierr /= SIONLIB_noerr) exit
         else if (trim(category) == 'cur') then
            call mpas_decode(handle % metadata_buffer, handle % metadata_cursor, test_start_cursor)
#ifdef SIONLIB_DEBUG
            call mpas_log_write('Reading start_cursor metadata for consistency checking ... ')
            call mpas_log_write('stored start_cursor      = $i', intArgs=(/ int(test_start_cursor) /) )
            call mpas_log_write('handle % start_cursor    = $i', intArgs=(/ int(handle % start_cursor) /) )
#endif
            if (test_start_cursor /= handle % start_cursor) then
               call mpas_log_write('ERROR while reading SIONlib metadata, stored and actual start cursor do not match: $i vs $i', &
                                   messageType=MPAS_LOG_ERR, intArgs=(/ int(test_start_cursor), int(handle % start_cursor) /) )
               ierr = SIONLIB_err
               exit
            end if
         else if (trim(category) == 'eof') then
#ifdef SIONLIB_DEBUG
            call mpas_log_write('Read eof statement, exit loop.')
#endif
            ierr = SIONLIB_noerr
            exit
#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
         else if (trim(category) == 'dbg') then
            prepend_check = .true.
#else
         else if (trim(category) == 'dbg') then
            call mpas_log_write('ERROR while reading SIONlib metadata, found prepend variable properites identifier "' // trim(category) // '" and prepending is switched off', messageType=MPAS_LOG_ERR)
            ierr = SIONLIB_err
            exit
#endif
         else
            call mpas_log_write('ERROR while reading SIONlib metadata, unknown identifier "' // trim(category) // '"', messageType=MPAS_LOG_ERR)
            ierr = SIONLIB_err
            exit
         end if
      end do

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      if (.not. prepend_check) then
         call mpas_log_write('ERROR while reading SIONlib metadata, prepend variable properties identifier "dbg" not found and prepending is switched on', messageType=MPAS_LOG_ERR)
         ierr = SIONLIB_err
      else
         call mpas_log_write('Required prepend variable properties identifier "dbg" found while reading SIONlib metadata')
      end if
#endif

   contains

      subroutine MPAS_SION_read_metadata_def_dim(handle, dimname, dimid, dimsize)

         implicit none

         type (MPAS_IO_Handle_type), intent(inout) :: handle
         character (len=*), intent(in) :: dimname
         integer, intent(in) :: dimid
         integer, intent(in) :: dimsize

         type (dimlist_type), pointer :: new_dimlist_node

         allocate(new_dimlist_node)
         nullify(new_dimlist_node % next)
         allocate(new_dimlist_node % dimhandle)

         new_dimlist_node % dimhandle % dimname = trim(dimname)
         new_dimlist_node % dimhandle % dimid   = dimid
         new_dimlist_node % dimhandle % dimsize = dimsize
         if (dimsize == MPAS_IO_UNLIMITED_DIM) then
            new_dimlist_node % dimhandle % is_unlimited_dim = .true.
            ! DH* check what happens with this in case
            ! netCDF is used directly instead of PIO *DH
            handle % unlimited_dimid = dimid
         else
            new_dimlist_node % dimhandle % is_unlimited_dim = .false.
         end if

         if (.not. associated(handle % dimlist_head)) then
            handle % dimlist_head => new_dimlist_node
            !call mpas_log_write('Assigning head for '//trim(dimname))
         end if
         if (.not. associated(handle % dimlist_tail)) then
            handle % dimlist_tail => new_dimlist_node
            !call mpas_log_write('Assigning tail for '//trim(dimname))
         else
            handle % dimlist_tail % next => new_dimlist_node
            handle % dimlist_tail => handle % dimlist_tail % next
            !call mpas_log_write('Extending tail for '//trim(dimname))
         end if

      end subroutine MPAS_SION_read_metadata_def_dim

      subroutine MPAS_SION_read_metadata_def_att(handle, attname, varid, atttype, intval, &
                                               & intarray, realval, realarray, charval, ierr)

         implicit none

         type (MPAS_IO_Handle_type), intent(inout) :: handle
         character (len=*), intent(in) :: attname
         integer, intent(in) :: varid
         integer, intent(in) :: atttype
         integer, intent(in), optional :: intval
         integer, dimension(:), intent(in), optional :: intarray
         real (KIND=RKIND), intent(in), optional :: realval
         real (KIND=RKIND), dimension(:), intent(in), optional :: realarray
         character (len=*), intent(in), optional :: charval
         integer, intent(out) :: ierr

         logical :: found
         type (attlist_type), pointer :: new_att_node

         allocate(new_att_node)
         nullify(new_att_node % next)
         allocate(new_att_node % atthandle)

         new_att_node % atthandle % attName = trim(attname)
         new_att_node % atthandle % fieldid = varid
         new_att_node % atthandle % attType = atttype

         if (present(intval)) then
            new_att_node % atthandle % attValueInt = intval
         else if (present(intarray)) then
            allocate(new_att_node % atthandle % attValueIntA(size(intarray)))
            new_att_node % atthandle % attValueIntA(:) = intarray(:)
         else if (present(realval)) then
            new_att_node % atthandle % attValueReal = realval
         else if (present(realarray)) then
            allocate(new_att_node % atthandle % attValueRealA(size(realarray)))
            new_att_node % atthandle % attValueRealA(:) = realarray(:)
         else if (present(charval)) then
            new_att_node % atthandle % attValueText = trim(charval)
         end if
         ! DH* - is this correct? I believe it is, see mpas_io_types.inc
         new_att_node % atthandle % precision = MPAS_IO_NATIVE_PRECISION

         if (varid == SIONLIB_global) then
            if (.not. associated(handle % attlist_head)) then
               handle % attlist_head => new_att_node
               !call mpas_log_write('Assigning handle -> att head for global attribute ' // trim(new_att_node %  atthandle % attname))
            end if
            if (.not. associated(handle % attlist_tail)) then
               handle % attlist_tail => new_att_node
               !call mpas_log_write('Assigning handle -> att tail for global attribute ' // trim(new_att_node %  atthandle % attname))
            else
               handle % attlist_tail % next => new_att_node
               handle % attlist_tail => handle % attlist_tail % next
               !call mpas_log_write('Extending handle -> att tail for global attribute ' // trim(new_att_node %  atthandle % attname))
            end if
         else
            found = .false.
            field_cursor => handle % fieldlist_head
            do while (associated(field_cursor))
               if (varid == field_cursor % fieldhandle % fieldid) then
                  found = .true.
                  exit
               end if
               field_cursor => field_cursor % next
            end do
            if (.not.found) then
               call mpas_log_write('ERROR, cannot assign attribute ' // trim(attname) // ' to non-existing variable with id $i', &
                              messageType=MPAS_LOG_ERR, intArgs=(/ varid /) )
               ierr = SIONLIB_err
               deallocate(new_att_node % atthandle)
               deallocate(new_att_node)
               return
            end if
            if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
               field_cursor % fieldhandle % attlist_head => new_att_node
               !call mpas_log_write('Assigning field_cursor % fieldhandle -> att head for variable attribute ' // trim(new_att_node % atthandle % attname))
            end if
            if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
               field_cursor % fieldhandle % attlist_tail => new_att_node
               !call mpas_log_write('Assigning field_cursor % fieldhandle -> att tail for variable attribute ' // trim(new_att_node % atthandle % attname))
            else
               field_cursor % fieldhandle % attlist_tail % next => new_att_node
               field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
               !call mpas_log_write('Extending field_cursor % fieldhandle -> att tail for variable attribute ' // trim(new_att_node % atthandle % attname))
            end if
         end if

         ierr = SIONLIB_noerr

      end subroutine MPAS_SION_read_metadata_def_att

      subroutine MPAS_SION_read_metadata_def_var(handle, varname, varid, vartype, ndims, dimids, ierr)

         implicit none

         type (MPAS_IO_Handle_type), intent(inout) :: handle
         character (len=*), intent(in) :: varname
         integer, intent(in) :: varid
         integer, intent(in) :: vartype
         integer, intent(in) :: ndims
         integer, dimension(:), intent(in) :: dimids
         integer, intent(out) :: ierr

         integer :: idim
         logical :: found
         type (fieldlist_type), pointer :: new_fieldlist_node

         allocate(new_fieldlist_node)
         nullify(new_fieldlist_node % next)
         allocate(new_fieldlist_node % fieldhandle)

         new_fieldlist_node % fieldhandle % fieldname  = trim(varname)
         new_fieldlist_node % fieldhandle % fieldid    = varid
         new_fieldlist_node % fieldhandle % field_type = vartype
         new_fieldlist_node % fieldhandle % ndims      = ndims

         if (vartype == MPAS_IO_DOUBLE) then
            new_fieldlist_node % fieldhandle % precision = MPAS_IO_DOUBLE_PRECISION
         else if (vartype == MPAS_IO_REAL) then
            new_fieldlist_node % fieldhandle % precision = MPAS_IO_SINGLE_PRECISION
         else
            new_fieldlist_node % fieldhandle % precision = MPAS_IO_NATIVE_PRECISION
         end if

         allocate(new_fieldlist_node % fieldhandle % dims(ndims))
         if (ndims > 0) then
            do idim=1,ndims
               new_fieldlist_node % fieldhandle % dims(idim) % dimid = dimids(idim)
               found = .false.
               dim_cursor => handle % dimlist_head
               do while (associated(dim_cursor))
                  if (dim_cursor % dimhandle % dimid == dimids(idim)) then
                     found = .true.
                     exit
                  end if
                  dim_cursor => dim_cursor % next
               end do
               if (.not.found) then
                  call mpas_log_write('ERROR, dimension with id $i, required for variable ' // trim(varname) // ' not found.', &
                                 messageType=MPAS_LOG_ERR, intArgs=(/ dimids(idim) /) )
                  ierr = SIONLIB_err
                  deallocate(new_fieldlist_node % fieldhandle)
                  deallocate(new_fieldlist_node)
                  return
               end if
               new_fieldlist_node % fieldhandle % dims(idim) % dimname = trim(dim_cursor % dimhandle % dimname)
               new_fieldlist_node % fieldhandle % dims(idim) % dimsize = dim_cursor % dimhandle % dimsize
               new_fieldlist_node % fieldhandle % dims(idim) % is_unlimited_dim = dim_cursor % dimhandle % is_unlimited_dim
               if (new_fieldlist_node % fieldhandle % dims(idim) % is_unlimited_dim) then
                  new_fieldlist_node % fieldhandle % has_unlimited_dim = .true.
               end if
            end do
         end if

         if (.not. associated(handle % fieldlist_head)) then
            handle % fieldlist_head => new_fieldlist_node
            !call mpas_log_write('Assigning head for '//trim(varname))
         end if
         if (.not. associated(handle % fieldlist_tail)) then
            handle % fieldlist_tail => new_fieldlist_node
            !call mpas_log_write('Assigning tail for '//trim(varname))
         else
            handle % fieldlist_tail % next => new_fieldlist_node
            handle % fieldlist_tail => handle % fieldlist_tail % next
            !call mpas_log_write('Extending tail for '//trim(varname))
         end if

         ierr = SIONLIB_noerr

      end subroutine MPAS_SION_read_metadata_def_var

      subroutine MPAS_SION_read_metadata_set_pos(handle, varid, position, varsize, chunksize, ierr)

         implicit none

         type (MPAS_IO_Handle_type), intent(inout) :: handle
         integer, intent(in) :: varid
         integer (kind=I8KIND), intent(in) :: position
         integer (kind=I8KIND), intent(in) :: varsize
         integer (kind=I8KIND), intent(in) :: chunksize
         integer, intent(out) :: ierr

         logical :: found
         type (fieldlist_type), pointer :: field_cursor

         found = .false.
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (varid == field_cursor % fieldhandle % fieldid) then
               found = .true.
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not.found) then
            call mpas_log_write('ERROR, cannot assign position to non-existing variable with id $i', &
                           messageType=MPAS_LOG_ERR, intArgs=(/ varid /) )
            ierr = SIONLIB_err
            return
         end if

         field_cursor % fieldhandle % sionlib_position  = position
         field_cursor % fieldhandle % sionlib_fieldsize = varsize
         field_cursor % fieldhandle % sionlib_chunksize = chunksize

         ierr = SIONLIB_noerr

      end subroutine MPAS_SION_read_metadata_set_pos

#endif
   end subroutine MPAS_SION_read_metadata

   subroutine MPAS_SION_read_array(handle, field_cursor, &
      & realval, realarray1d, realarray2d, realarray3d, realarray4d, realarray5d, &
      & intval, intarray1d, intarray2d, intarray3d, intarray4d, intarray5d, &
      & singleval, singlearray1d, singlearray2d, singlearray3d, singlearray4d, singlearray5d, &
      & doubleval, doublearray1d, doublearray2d, doublearray3d, doublearray4d, doublearray5d, &
      & charval, chararray1d, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      type (fieldlist_type), pointer, intent(in) :: field_cursor
      real (kind=RKIND),  intent(out), optional :: realval
      real (kind=RKIND),  dimension(:), intent(out), optional :: realarray1d
      real (kind=RKIND),  dimension(:,:), intent(out), optional :: realarray2d
      real (kind=RKIND),  dimension(:,:,:), intent(out), optional :: realarray3d
      real (kind=RKIND),  dimension(:,:,:,:), intent(out), optional :: realarray4d
      real (kind=RKIND),  dimension(:,:,:,:,:), intent(out), optional :: realarray5d
      integer,            intent(out), optional :: intval
      integer,            dimension(:), intent(out), optional :: intarray1d
      integer,            dimension(:,:), intent(out), optional :: intarray2d
      integer,            dimension(:,:,:), intent(out), optional :: intarray3d
      integer,            dimension(:,:,:,:), intent(out), optional :: intarray4d
      integer,            dimension(:,:,:,:,:), intent(out), optional :: intarray5d
      real (kind=R4KIND), intent(out), optional :: singleval
      real (kind=R4KIND), dimension(:), intent(out), optional :: singlearray1d
      real (kind=R4KIND), dimension(:,:), intent(out), optional :: singlearray2d
      real (kind=R4KIND), dimension(:,:,:), intent(out), optional :: singlearray3d
      real (kind=R4KIND), dimension(:,:,:,:), intent(out), optional :: singlearray4d
      real (kind=R4KIND), dimension(:,:,:,:,:), intent(out), optional :: singlearray5d
      real (kind=R8KIND), intent(out), optional :: doubleval
      real (kind=R8KIND), dimension(:), intent(out), optional :: doublearray1d
      real (kind=R8KIND), dimension(:,:), intent(out), optional :: doublearray2d
      real (kind=R8KIND), dimension(:,:,:), intent(out), optional :: doublearray3d
      real (kind=R8KIND), dimension(:,:,:,:), intent(out), optional :: doublearray4d
      real (kind=R8KIND), dimension(:,:,:,:,:), intent(out), optional :: doublearray5d
      character (len=*),  intent(out), optional :: charval
      character (len=*),  dimension(:), intent(out), optional :: chararray1d
      integer, intent(out) :: ierr

#ifdef SIONLIB
      integer   :: varid
      integer   :: vartype
      integer (kind=I8KIND) :: position
      integer (kind=I8KIND) :: varsize
      integer (kind=I8KIND) :: chunksize
      integer :: checksum
      integer (kind=I8KIND) :: count
      character (len=StrKIND) :: padded_charval
      character (len=HugeStrKIND) :: padded_charval_huge
      integer :: ichar

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      integer   :: varid_test
      integer   :: vartype_test
      integer (kind=I8KIND) :: varsize_test
#endif

      varid     = field_cursor % fieldhandle % fieldid
      vartype   = field_cursor % fieldhandle % field_type
      position  = field_cursor % fieldhandle % sionlib_position
      varsize   = field_cursor % fieldhandle % sionlib_fieldsize
      chunksize = field_cursor % fieldhandle % sionlib_chunksize

      ! Calculated and stored positions are without considering the start cursor position of each task
      position = position + handle % start_cursor
      call MPAS_SION_set_cursor(handle, position, ierr)
      if (ierr /= SIONLIB_noerr) return

#ifdef SIONLIB_PREPEND_VARIABLE_PROPERTIES
      call FSION_COLL_FREAD_MPI(varid_test,   int(kind(varid_test),8),   int(1,8), handle % sidd, count)
      call FSION_COLL_FREAD_MPI(vartype_test, int(kind(vartype_test),8), int(1,8), handle % sidd, count)
      call FSION_COLL_FREAD_MPI(varsize_test, int(kind(varsize_test),8), int(1,8), handle % sidd, count)
      if (varid_test /= varid) then
         call mpas_log_write('ERROR mismatch of varid, expected $i vs. read $i.', messageType=MPAS_LOG_ERR, intArgs=(/ varid, varid_test /) )
         ierr = SIONLIB_err
         return
      else if (vartype_test /= vartype) then
         call mpas_log_write('ERROR mismatch of vartype, expected $i vs. read $i.', messageType=MPAS_LOG_ERR, intArgs=(/ vartype, vartype_test /) )
         ierr = SIONLIB_err
         return
      else if (varsize_test /= varsize) then
         call mpas_log_write('ERROR mismatch of varsize, expected $i vs. read $i.', messageType=MPAS_LOG_ERR, intArgs=(/ int(varsize), int(varsize_test) /) )
         ierr = SIONLIB_err
         return
      end if
#endif

      if (present(realval)) then
         checksum = 1
         call FSION_COLL_FREAD_MPI(realval, int(kind(realval),8), int(1,8), handle % sidd, count)
      else if (present(realarray1d)) then
         checksum = size(realarray1d)
         call FSION_COLL_FREAD_MPI(realarray1d(1), int(kind(realarray1d(1)),8), int(size(realarray1d),8), handle % sidd, count)
      else if (present(realarray2d)) then
         checksum = size(realarray2d)
         call FSION_COLL_FREAD_MPI(realarray2d(1,1), int(kind(realarray2d(1,1)),8), int(size(realarray2d),8), handle % sidd, count)
      else if (present(realarray3d)) then
         checksum = size(realarray3d)
         call FSION_COLL_FREAD_MPI(realarray3d(1,1,1), int(kind(realarray3d(1,1,1)),8), int(size(realarray3d),8), handle % sidd, count)
      else if (present(realarray4d)) then
         checksum = size(realarray4d)
         call FSION_COLL_FREAD_MPI(realarray4d(1,1,1,1), int(kind(realarray4d(1,1,1,1)),8), int(size(realarray4d),8), handle % sidd, count)
      else if (present(realarray5d)) then
         checksum = size(realarray5d)
         call FSION_COLL_FREAD_MPI(realarray5d(1,1,1,1,1), int(kind(realarray5d(1,1,1,1,1)),8), int(size(realarray5d),8), handle % sidd, count)
      else if (present(intval)) then
         checksum = 1
         call FSION_COLL_FREAD_MPI(intval, int(kind(intval),8), int(1,8), handle % sidd, count)
      else if (present(intarray1d)) then
         checksum = size(intarray1d)
         call FSION_COLL_FREAD_MPI(intarray1d(1), int(kind(intarray1d(1)),8), int(size(intarray1d),8), handle % sidd, count)
      else if (present(intarray2d)) then
         checksum = size(intarray2d)
         call FSION_COLL_FREAD_MPI(intarray2d(1,1), int(kind(intarray2d(1,1)),8), int(size(intarray2d),8), handle % sidd, count)
      else if (present(intarray3d)) then
         checksum = size(intarray3d)
         call FSION_COLL_FREAD_MPI(intarray3d(1,1,1), int(kind(intarray3d(1,1,1)),8), int(size(intarray3d),8), handle % sidd, count)
      else if (present(intarray4d)) then
         checksum = size(intarray4d)
         call FSION_COLL_FREAD_MPI(intarray4d(1,1,1,1), int(kind(intarray4d(1,1,1,1)),8), int(size(intarray4d),8), handle % sidd, count)
      else if (present(intarray5d)) then
         checksum = size(intarray5d)
         call FSION_COLL_FREAD_MPI(intarray5d(1,1,1,1,1), int(kind(intarray5d(1,1,1,1,1)),8), int(size(intarray5d),8), handle % sidd, count)
      else if (present(singleval)) then
         checksum = 1
         call FSION_COLL_FREAD_MPI(singleval, int(kind(singleval),8), int(1,8), handle % sidd, count)
      else if (present(singlearray1d)) then
         checksum = size(singlearray1d)
         call FSION_COLL_FREAD_MPI(singlearray1d(1), int(kind(singlearray1d(1)),8), int(size(singlearray1d),8), handle % sidd, count)
      else if (present(singlearray2d)) then
         checksum = size(singlearray2d)
         call FSION_COLL_FREAD_MPI(singlearray2d(1,1), int(kind(singlearray2d(1,1)),8), int(size(singlearray2d),8), handle % sidd, count)
      else if (present(singlearray3d)) then
         checksum = size(singlearray3d)
         call FSION_COLL_FREAD_MPI(singlearray3d(1,1,1), int(kind(singlearray3d(1,1,1)),8), int(size(singlearray3d),8), handle % sidd, count)
      else if (present(singlearray4d)) then
         checksum = size(singlearray4d)
         call FSION_COLL_FREAD_MPI(singlearray4d(1,1,1,1), int(kind(singlearray4d(1,1,1,1)),8), int(size(singlearray4d),8), handle % sidd, count)
      else if (present(singlearray5d)) then
         checksum = size(singlearray5d)
         call FSION_COLL_FREAD_MPI(singlearray5d(1,1,1,1,1), int(kind(singlearray5d(1,1,1,1,1)),8), int(size(singlearray5d),8), handle % sidd, count)
      else if (present(doubleval)) then
         checksum = 1
         call FSION_COLL_FREAD_MPI(doubleval, int(kind(doubleval),8), int(1,8), handle % sidd, count)
      else if (present(doublearray1d)) then
         checksum = size(doublearray1d)
         call FSION_COLL_FREAD_MPI(doublearray1d(1), int(kind(doublearray1d(1)),8), int(size(doublearray1d),8), handle % sidd, count)
      else if (present(doublearray2d)) then
         checksum = size(doublearray2d)
         call FSION_COLL_FREAD_MPI(doublearray2d(1,1), int(kind(doublearray2d(1,1)),8), int(size(doublearray2d),8), handle % sidd, count)
      else if (present(doublearray3d)) then
         checksum = size(doublearray3d)
         call FSION_COLL_FREAD_MPI(doublearray3d(1,1,1), int(kind(doublearray3d(1,1,1)),8), int(size(doublearray3d),8), handle % sidd, count)
      else if (present(doublearray4d)) then
         checksum = size(doublearray4d)
         call FSION_COLL_FREAD_MPI(doublearray4d(1,1,1,1), int(kind(doublearray4d(1,1,1,1)),8), int(size(doublearray4d),8), handle % sidd, count)
      else if (present(doublearray5d)) then
         checksum = size(doublearray5d)
         call FSION_COLL_FREAD_MPI(doublearray5d(1,1,1,1,1), int(kind(doublearray5d(1,1,1,1,1)),8), int(size(doublearray5d),8), handle % sidd, count)
      else if (present(charval)) then
         checksum = len(padded_charval)
         call FSION_COLL_FREAD_MPI(padded_charval, int(1,8), int(len(padded_charval),8), handle % sidd, count)
         charval = trim(padded_charval)
      else if (present(chararray1d)) then
         checksum = size(chararray1d)*len(padded_charval)
         if (checksum > len(padded_charval_huge)) then
             call mpas_log_write('MPAS_SION_read_darray: error, buffer overflow for reading chararray1d from file for variable $i - increase HugeStrKIND!', &
                                 messageType=MPAS_LOG_ERR, intArgs=(/ field_cursor % fieldhandle % fieldid /) )
             ierr = SIONLIB_err
             return
         end if
         ! Read entire array as one long string and split up
         call FSION_COLL_FREAD_MPI(padded_charval_huge, int(1,8), int(len(padded_charval_huge),8), handle % sidd, count)
         do ichar=1,size(chararray1d)
            padded_charval = padded_charval_huge((ichar-1)*len(padded_charval)+1:ichar*len(padded_charval))
            chararray1d(ichar) = trim(padded_charval)
         end do
         ! DH*
         call mpas_log_write('Attention, reading charArray1d variables has not been tested yet!', messageType=MPAS_LOG_WARN)
         ! *DH
      else
          call mpas_log_write('MPAS_SION_read_darray: no known storage type specified for variable $i', &
                              messageType=MPAS_LOG_ERR, intArgs=(/ varid /) )
          ierr = SIONLIB_err
          return
      end if

      if (varsize /= checksum .or. count /= checksum) then
         call mpas_log_write('Error in MPAS_SION_read_array - size mismatch while reading data for varid $i: arsize, count, checksum are $i, $i, $i', &
                             messageType=MPAS_LOG_ERR, intArgs=(/ varid, int(varsize), int(count), checksum /) )
         ierr = SIONLIB_err
         return
      end if

      ierr = SIONLIB_noerr
#endif

   end subroutine MPAS_SION_read_array

   subroutine MPAS_SION_set_cursor(handle, position, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer (kind=I8KIND), intent(in) :: position
      integer, intent(out) :: ierr

#ifdef SIONLIB
      integer (kind=I8KIND) :: current

      current = FSION_GET_POSITION(handle % sidd)
      if (current == position) then
         ! Already at correct location, nothing to do
         ierr = SIONLIB_noerr
         return
      end if
#ifdef SIONLIB_DEBUG
      call mpas_log_write('Positioning cursor at target $i; current position and start block are $i and $i', &
                          intArgs=(/ int(position), int(current), int(handle % start_cursor) /) )
#endif

      call FSION_SEEK(handle % sidd, SION_CURRENT_RANK, SION_CURRENT_BLK, (position - handle % start_cursor), ierr)
      if (ierr/=1) then
         call mpas_log_write('ERROR in MPAS_SION_set_cursor - FSION_SEEK failed', messageType=MPAS_LOG_ERR)
         ierr = SIONLIB_err
         return
      end if
      current = FSION_GET_POSITION(handle % sidd)
      if (current /= position) then
         call mpas_log_write('ERROR in MPAS_SION_set_cursor - after seek, current position != target position: $i != $i', &
                             messageType=MPAS_LOG_ERR, intArgs=(/ int(current), int(position) /) )
         ierr = SIONLIB_err
         return
      end if
      ierr = SIONLIB_noerr
#endif

      !DIR$ ATTRIBUTES INLINE :: MPAS_SION_set_cursor
   end subroutine MPAS_SION_set_cursor

end module mpas_sion
