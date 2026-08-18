!     *****************************************************************
!     **                                                             **
!     **                 CKINTP capacity preflight                   **
!     **                                                             **
!     **  Reject mechanisms that exceed the fixed legacy CKINTP      **
!     **  work-array/string limits before CKINTP mutates parser      **
!     **  state or writes a linking file.                            **
!     **                                                             **
!     *****************************************************************

      subroutine ckintp_preflight

      use chemistry_setup, only: mechdir
      use chemistry_string_limits, only: species_name_len, &
                                          mechanism_line_len
      use iso_fortran_env, only: output_unit

      implicit none

      integer, parameter :: ckintp_species_capacity = 3000
      integer, parameter :: ckintp_reaction_capacity = 9000
      integer, parameter :: ckintp_max_participants = 8
      integer, parameter :: scan_line_len = 4096
      integer, parameter :: max_tokens = 128
      integer, parameter :: section_none = 0
      integer, parameter :: section_species = 1
      integer, parameter :: section_reactions = 2

      character(len=scan_line_len) :: raw, line, equation
      character(len=mechanism_line_len) :: tokens(max_tokens)
      character(len=64) :: keyword
      character(len=512) :: mechfile, linkfile, outfile
      integer :: lun, ios, line_number, section
      integer :: species_count, reaction_count
      integer :: max_name_len, max_participants, max_record_len
      integer :: first_long_name_line, first_many_participants_line
      integer :: first_long_record_line, first_reaction_overflow_line
      integer :: ntokens, i, participants, lout
      logical :: is_reaction, failed

      mechfile = trim(mechdir)//'chem.inp'
      linkfile = trim(mechdir)//'cklink'
      outfile = trim(mechdir)//'chem.out'

      species_count = 0
      reaction_count = 0
      max_name_len = 0
      max_participants = 0
      max_record_len = 0
      first_long_name_line = 0
      first_many_participants_line = 0
      first_long_record_line = 0
      first_reaction_overflow_line = 0
      line_number = 0
      section = section_none

      open(newunit=lun, file=trim(mechfile), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) return

      do
         read(lun, '(A)', iostat=ios) raw
         if (ios < 0) exit
         if (ios > 0) then
            close(lun)
            write(output_unit, '(A,I0,A,A)') &
               ' Error...CKINTP preflight could not read line ', &
               line_number + 1, ' of ', trim(mechfile)
            call remove_file_if_present(linkfile)
            error stop 1
         endif

         line_number = line_number + 1
         call strip_comment(raw, line)
         if (len_trim(line) == 0) cycle

         call tokenize(line, tokens, ntokens)
         if (ntokens == 0) cycle
         keyword = upper(tokens(1))

         select case (trim(keyword))
         case ('SPECIES', 'SPEC')
            section = section_species
            if (ntokens > 1) then
               do i = 2, ntokens
                  call record_species(tokens(i), line_number)
               enddo
            endif
            cycle
         case ('REACTIONS', 'REAC')
            section = section_reactions
            cycle
         case ('ELEMENTS', 'ELEM', 'THERMO', 'THER')
            section = section_none
            cycle
         case ('END')
            if (section == section_reactions) exit
            section = section_none
            cycle
         end select

         select case (section)
         case (section_species)
            do i = 1, ntokens
               call record_species(tokens(i), line_number)
            enddo

         case (section_reactions)
            max_record_len = max(max_record_len, len_trim(line))
            if (len_trim(line) > mechanism_line_len .and. &
                first_long_record_line == 0) then
               first_long_record_line = line_number
            endif

            is_reaction = index(line, '/') == 0 .and. &
                          index(upper(tokens(1)), 'DUP') /= 1
            if (.not. is_reaction) cycle

            reaction_count = reaction_count + 1
            if (reaction_count == ckintp_reaction_capacity + 1) then
               first_reaction_overflow_line = line_number
            endif

            call extract_equation(tokens, ntokens, equation)
            participants = reaction_participants(equation)
            max_participants = max(max_participants, participants)
            if (participants > ckintp_max_participants .and. &
                first_many_participants_line == 0) then
               first_many_participants_line = line_number
            endif
         end select
      enddo
      close(lun)

      failed = species_count > ckintp_species_capacity .or. &
               reaction_count > ckintp_reaction_capacity .or. &
               max_name_len > species_name_len .or. &
               max_participants > ckintp_max_participants .or. &
               max_record_len > mechanism_line_len

      if (.not. failed) return

      call remove_file_if_present(linkfile)
      call write_report(output_unit)
      flush(output_unit)

      open(newunit=lout, file=trim(outfile), status='replace', &
           action='write', iostat=ios)
      if (ios == 0) then
         call write_report(lout)
         close(lout)
      endif

      error stop 1

      contains

      subroutine write_report(unit_number)
      integer, intent(in) :: unit_number

      write(unit_number, '(A)') &
         ' Error...CKINTP preflight rejected the mechanism before parsing.'
      if (reaction_count > ckintp_reaction_capacity) then
         write(unit_number, '(A,I0,A,I0,A,I0,A)') &
            '   reactions: ', reaction_count, ' (capacity IDIM=', &
            ckintp_reaction_capacity, &
            '; first overflow at input line ', &
            first_reaction_overflow_line, ')'
      endif
      if (species_count > ckintp_species_capacity) then
         write(unit_number, '(A,I0,A,I0,A)') &
            '   species: ', species_count, ' (capacity KDIM=', &
            ckintp_species_capacity, ')'
      endif
      if (max_participants > ckintp_max_participants) then
         write(unit_number, '(A,I0,A,I0,A,I0,A)') &
            '   max reaction participants: ', max_participants, &
            ' (capacity MAXSP=', ckintp_max_participants, &
            '; first overflow at input line ', &
            first_many_participants_line, ')'
      endif
      if (max_name_len > species_name_len) then
         write(unit_number, '(A,I0,A,I0,A,I0,A)') &
            '   max species-name length: ', max_name_len, &
            ' (capacity LSYM=', species_name_len, &
            '; first overflow at input line ', first_long_name_line, ')'
      endif
      if (max_record_len > mechanism_line_len) then
         write(unit_number, '(A,I0,A,I0,A,I0,A)') &
            '   max mechanism-record length: ', max_record_len, &
            ' (capacity=', mechanism_line_len, &
            '; first overflow at input line ', &
            first_long_record_line, ')'
      endif
      write(unit_number, '(A)') &
         '   No reaction or auxiliary records were passed to CKINTP.'
      write(unit_number, '(A)') &
         '   No cklink file was produced.'
      end subroutine write_report


      subroutine strip_comment(input, output)
      character(len=*), intent(in) :: input
      character(len=*), intent(out) :: output
      integer :: bang

      output = input
      bang = index(output, '!')
      if (bang > 0) output(bang:) = ' '
      output = adjustl(output)
      end subroutine strip_comment


      subroutine tokenize(input, output_tokens, output_count)
      character(len=*), intent(in) :: input
      character(len=*), intent(out) :: output_tokens(:)
      integer, intent(out) :: output_count
      integer :: first, last, n

      output_tokens = ' '
      output_count = 0
      n = len_trim(input)
      first = 1
      do while (first <= n)
         do while (first <= n .and. is_space(input(first:first)))
            first = first + 1
         enddo
         if (first > n) exit
         last = first
         do while (last <= n .and. .not. is_space(input(last:last)))
            last = last + 1
         enddo
         if (output_count == size(output_tokens)) exit
         output_count = output_count + 1
         output_tokens(output_count) = input(first:last-1)
         first = last + 1
      enddo
      end subroutine tokenize


      subroutine record_species(name, input_line)
      character(len=*), intent(in) :: name
      integer, intent(in) :: input_line
      integer :: n

      n = len_trim(name)
      if (n == 0) return
      species_count = species_count + 1
      max_name_len = max(max_name_len, n)
      if (n > species_name_len .and. first_long_name_line == 0) then
         first_long_name_line = input_line
      endif
      end subroutine record_species


      subroutine extract_equation(input_tokens, input_count, output)
      character(len=*), intent(in) :: input_tokens(:)
      integer, intent(in) :: input_count
      character(len=*), intent(out) :: output
      integer :: item, equation_tokens, parse_status
      double precision :: dummy

      output = ' '
      equation_tokens = input_count
      if (input_count >= 4) then
         read(input_tokens(input_count), *, iostat=parse_status) dummy
         if (parse_status == 0) then
            read(input_tokens(input_count-1), *, &
                 iostat=parse_status) dummy
            if (parse_status == 0) then
               read(input_tokens(input_count-2), *, &
                    iostat=parse_status) dummy
               if (parse_status == 0) equation_tokens = input_count - 3
            endif
         endif
      endif

      do item = 1, equation_tokens
         if (len_trim(output) + len_trim(input_tokens(item)) > &
             len(output)) exit
         output = trim(output)//trim(input_tokens(item))
      enddo
      end subroutine extract_equation


      integer function reaction_participants(input) result(total)
      character(len=*), intent(in) :: input
      integer :: arrow_pos, arrow_len

      total = 0
      call find_arrow(input, arrow_pos, arrow_len)
      if (arrow_pos == 0) return
      total = count_side_terms(input(:arrow_pos-1)) + &
              count_side_terms(input(arrow_pos+arrow_len:))
      end function reaction_participants


      subroutine find_arrow(input, arrow_pos, arrow_len)
      character(len=*), intent(in) :: input
      integer, intent(out) :: arrow_pos, arrow_len

      arrow_pos = index(input, '<=>')
      if (arrow_pos > 0) then
         arrow_len = 3
         return
      endif
      arrow_pos = index(input, '=>')
      if (arrow_pos > 0) then
         arrow_len = 2
         return
      endif
      arrow_pos = index(input, '=')
      if (arrow_pos > 0) then
         arrow_len = 1
      else
         arrow_len = 0
      endif
      end subroutine find_arrow


      integer function count_side_terms(side) result(count)
      character(len=*), intent(in) :: side
      character(len=len(side)) :: term
      integer :: item, first, depth, n

      count = 0
      n = len_trim(side)
      if (n == 0) return
      first = 1
      depth = 0
      do item = 1, n + 1
         if (item <= n) then
            select case (side(item:item))
            case ('(')
               depth = depth + 1
            case (')')
               depth = max(0, depth - 1)
            case ('+')
               if (depth /= 0) cycle
               term = ' '
               if (item > first) term = side(first:item-1)
               if (.not. is_third_body(term)) count = count + 1
               first = item + 1
            end select
         else
            term = ' '
            if (n >= first) term = side(first:n)
            if (.not. is_third_body(term)) count = count + 1
         endif
      enddo
      end function count_side_terms


      logical function is_third_body(term) result(is_tb)
      character(len=*), intent(in) :: term
      character(len=len(term)) :: cleaned

      cleaned = upper(adjustl(term))
      is_tb = trim(cleaned) == 'M' .or. trim(cleaned) == 'HV' .or. &
              len_trim(cleaned) == 0
      end function is_third_body


      logical function is_space(ch)
      character(len=1), intent(in) :: ch

      is_space = ch == ' ' .or. iachar(ch) == 9
      end function is_space


      pure function upper(text) result(output)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: output
      integer :: item, code

      output = text
      do item = 1, len(text)
         code = iachar(output(item:item))
         if (code >= iachar('a') .and. code <= iachar('z')) then
            output(item:item) = &
               achar(code - iachar('a') + iachar('A'))
         endif
      enddo
      end function upper


      subroutine remove_file_if_present(path)
      character(len=*), intent(in) :: path
      integer :: delete_unit, delete_status
      logical :: exists

      inquire(file=trim(path), exist=exists)
      if (.not. exists) return
      open(newunit=delete_unit, file=trim(path), status='old', &
           action='readwrite', iostat=delete_status)
      if (delete_status == 0) close(delete_unit, status='delete')
      end subroutine remove_file_if_present

      end subroutine ckintp_preflight
