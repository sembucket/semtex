% [0.1j]   doubled commas after title and superfluous comma before brackets 
%                removed (resulted by last changes), reported by Maarten Wisse
% [0.1i]   some warnings for missing authors/editors corrected, reported by Maarten Wisse
% [0.1h]   changed handling of authors/editors in @book, reported by Maaarten Wisse:
%		The golden rule with regard to editors and authors is: 
%			1. no author but only editor: then the 
%	   		   editor goes before the title. 
%			2. if there is an author \textbf{and} 
%			   an editor, then the \emph{author} goes
%			   \textbf{before} the title, and the 
%			   \emph{editor} goes \textbf{after} it.
% [0.1g]   superfluous \bibatsep in `format.title' removed, reported by Maarten Wisse
% [0.1f]   forgotten argument of \jbbibargs inserted
% [0.1e]   \bibbudcsep introduced instead of `~-- '
% [0.1d]   \bibedinformat introduced for inline editors
% [0.1c]   \el for bibliography renamed to \bibel, some stuff for identifying authors added
% [0.1b]   format.pages added for @BOOK and @BOOKLET
%          for names which should not be recognized by BibTex, use 
%               `author = {Gregory{ }the{ }Great}'
%           instead of
%               `author = {{Gregory the Great}}'
%          because of special functions for capitalizing the `von' part.
% [0.1a]   \bibAnnote stuff added
% [0.1]    jox.bst has reached beta stage 
%
% [0.0.1p] spaces after control sequences added for better linebreaking 
%            regarding the well known bibtex bug
% [0.0.1o] \volname added
% [0.0.1n] @incollection: `booktitleaddon' follows booktitle and volume, preceds series;
%          no \volumename if `volume' comes with citation of an essay;
%          `capitalize.von.part' added (Thanks to Oren Patashnik)
% [0.0.1m] `booktitleaddon' follows immediately booktitle, \volumename inserted in 
%           `format.in.editor.booktitle.volume.voltitle'
% [0.0.1l] bug in `format.in.editor.booktitle.volume.voltitle' removed (completely rewritten, 
%           error message, if no volumetitle was given); a \bpubaddr replaced by \bibbdsep in `format.publisher';
%           `format.titleaddon' added for @unpublished
% [0.0.1k] experimental noed field added, noed=1 disables output of `(ed.)', removing `textedition'
% [0.0.1j] bug in `format.in.editor.booktitle.volume.voltitle' removed (error message, if no editor was given)
% [0.0.1i] `format.endnote' added to make the `endnote' option work
% [0.0.1h] `format.in.editor.booktitle.volume.voltitle' introduced for non crossref'd incollections, 
%           this makes `volformat=reversed' working for booktitle too
% [0.0.1g] `format.bvolume.btitle' for non crossref'd incollections added
% [0.0.1f] `volume' removed if volumetitle missing for incollection and inbook crossrefs
% [0.0.1e] missing year (if publisher is not given) in format.publisher.* added, volformat=reversed implemented, 
%          chapter for inbook moved back to pages 
% [0.0.1d] incollection adapted to give the same output with or without crossref
% [0.0.1c] comma before url for all types except @www
% [0.0.1b] First try to adapt the @in* types
% [0.0.1a] \artvolnumformat introduced, all types except the @in* adapted (oaddress.opublisher.oyear, make.translator, pseudonym etc.)

FUNCTION {version} { " v0.1j" }

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    stitle
    title
    type
    volume
    year
% 
    url
    urldate
    volumetitle
    juratitle
    shorttitle
    juraauthor
    shortauthor
    howcited
    ssedition
    language
    titleaddon
    booktitleaddon
    editortype
    sortkey
    annote
    gender
    translator
    flanguage
    oyear
    oaddress
    opublisher
    pseudonym
    volformat
    noed
%        
  }
  { repeat-authors }
%   ^^^^^^^^^^^^^^ JB  
  { label dummy.label extra.label sort.label short.list current-authors }
%         ^^^^^^^^^^^^^^^^^^^^^^^                       ^^^^^^^^^^^^^^^JB

INTEGERS { output.state before.all mid.sentence }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
}

STRINGS { s t u v w x }
%             ^^^^^^^ JB 

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }    
    { write$ mid.sentence 'output.state := }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
     { pop$ "empty " t * " in " * cite$ * warning$ }
     'output.nonnull
  if$
}

FUNCTION {format.endnote}
{ note empty$
      { "" }
      { "  \jbendnote {" note * "}" * }
  if$
}

%% added JB
%  FUNCTION {format.howcited}{ before.all 'output.state := "} {\bibhowcited}"  * }
FUNCTION {format.howcited}{ before.all 'output.state := "} {\bibhowcited}" format.endnote * * }

%% added JB, but not needed in oxford
% FUNCTION {format.howcited.artperiod}
% { before.all 'output.state := "} {\bibartperiodhowcited}" }

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.name.sig}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
  { s nameptr "{ff}{ vv}{ ll}{ jj}" format.name$ 't :=
      nameptr #1 >
        { namesleft #1 >
            { "" * t * }
            {  t "others" =
                { " et al." * }
                { "" * t * }
              if$
            }
          if$
        }
        't
      if$
    nameptr #1 + 'nameptr :=
    namesleft #1 -  'namesleft :=
  }
  while$
}

FUNCTION {format.snd.editors}
{ editor empty$
   { "" * }
   { author empty$
        { editortype empty$
             { editor num.names$ #1 >
                 { "\snded {\sndeditorsname}" * }
                 { "\snded {\sndeditorname}" * }
               if$
             }
             { "\snded {(" editortype * ")}" * * }
          if$
        }
        'skip$
    if$
   }
  if$
}

%% added JB
FUNCTION {set.after.author.sep}
{ before.all 'output.state := 
  author missing$
     { editor missing$
          { organization missing$
                 { "} {empty} {} {" }
% REMOVED                          ^ "*"
                 { "} {" author format.name.sig * "} {au} {" * }
            if$
          }
          { "} {" editor format.name.sig * "} {ed} {" * }
      if$
     }
     { "} {" author format.name.sig * "} {au} {" * }
  if$
  author empty$
      { format.snd.editors }
      'skip$
  if$
  "" *  
%  "} {" "" * 
}

%% added JB
FUNCTION {fin.entry}
{ format.howcited
  " \jbdoitem" * write$
  annote missing$
    { " \bibAnnoteFile {" cite$ * "}" * write$ newline$ }
    { " \bibAnnote {" annote * "}" * write$ newline$ }
  if$
}

%% added JB
FUNCTION {format.note}
{ note empty$
      'skip$
      { output.state before.all =
          { "\jbnote {}{" * note * "}" * }
          { " \jbnote {1}{" * note * "}" * }
       if$
      }
  if$
}

FUNCTION {add.blank}{  " " * before.all 'output.state := }

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {non.stop}
{ duplicate$
   "}" * %add.period$
   #-1 #1 substring$ "." =
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {capitalize}
{ "u" change.case$ "t" change.case$ }

FUNCTION {space.word}
{ " " swap$ * " " * }

FUNCTION {bbl.edition}
{ "\editionname {} " }

FUNCTION {bbl.vol}
{ "\volname {}" }

FUNCTION {bbl.volume}
{ "\volumename {}" }

FUNCTION {bbl.of}
{ " \volumeofname {} " }

FUNCTION {bbl.number}
{ "\numbername" }

FUNCTION {bbl.in}
{ "\incollinname {} " }

FUNCTION {bbl.in.series}
{ "\inseriesname {}" }

FUNCTION {bbl.chapter}
{ "\bibchaptername {}" }

FUNCTION {bbl.techrep}
{ "\technicalreportname {}" }

FUNCTION {bbl.mthesis}
{ "\mastersthesisname {}" }

FUNCTION {bbl.phdthesis}
{ "\phdthesisname {}"}

% added + changed JB removed "{} " to make string solutions possible
MACRO {jan} {"\janname"}
MACRO {feb} {"\febname"}
MACRO {mar} {"\marname"}
MACRO {apr} {"\aprname"}
MACRO {may} {"\mayname"}
MACRO {jun} {"\junname"}
MACRO {jul} {"\julname"}
MACRO {aug} {"\augname"}
MACRO {sep} {"\sepname"}
MACRO {oct} {"\octname"}
MACRO {nov} {"\novname"}
MACRO {dec} {"\decname"}

% added JB
FUNCTION {checked.emph.format}
{ type$ "book" = type$ "inbook" = type$ "proceedings" = type$ "inproceedings" = or or or
      { "\bibtfont {" }
      { "\bibapifont {"  }
  if$
}

% added JB
FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { checked.emph.format swap$ * "}" * 
%      type$ "book" = { "\bhook" * } 'skip$ if$ 
    }
  if$
}

% added JB, for booktitle in incollection and inproceedings if crossref's are missing
FUNCTION {b.emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\bibbtfont {" swap$ * "}" * }
  if$
}

% added JB
FUNCTION {write.url}
{ url empty$
      { urldate empty$
            'skip$
            { "there's an urldate, but no url in " cite$ * warning$ }               % warning added JB
        if$
      }
      { type$ "www" =
           'skip$
           { ", " * }
        if$
        title empty$ output.state before.all = and                                  % added JB to avoid leading whitespace for minimal misc and www
            { "" * mid.sentence 'output.state := }                                  % added JB to avoid leading whitespace for minimal misc and www
            { " " * }                                                               % added JB to avoid leading whitespace for minimal misc and www
        if$                                                                         % added JB to avoid leading whitespace for minimal misc and www
% comma before URL??????
        urldate empty$
                { "\jburluse {" cite$ * "}" * * }
%                { "\jburluse {" cite$ * "}~-- \urldatecomment " * urldate * * }
                { "\jburluse {" cite$ * "}\bibbudcsep {}\urldatecomment " * urldate * * }
        if$
      }
  if$
}

FUNCTION {lastnamemarkup}
{ duplicate$ empty$
     { pop$ "" }
     { "\bibnf {" swap$ * "}" * }
  if$
}

FUNCTION {editorlastnamemarkup}
{ duplicate$ empty$
     { pop$ "" }
     { author empty$
         { "\biblenf {" swap$ * "}" * } % lonely editor(s)!
         { "\bibenf {" swap$ * "}" * }
       if$
     }
  if$
}

FUNCTION {othermarkup}
{ duplicate$ empty$
       { pop$ " {}" }
       { " {" swap$ * "}" * }
  if$
}

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
  { s nameptr "{ll}" format.name$ 't := t lastnamemarkup 't :=
    s nameptr "{ff}" format.name$ othermarkup 'u :=
    s nameptr "{f.}" format.name$ othermarkup 'v :=
    s nameptr "{vv}" format.name$ othermarkup 'w :=
    s nameptr "{jj}" format.name$ othermarkup 'x :=
    t u * v * w * x * 't :=
     nameptr #1 >
      { namesleft #1 >
          { numnames #3 >
                { "" * }
                { numnames #3 <
                     { "\Bibbtasep " * t * }      % between (exactly) two authors
                     { "\Bibbfsasep " * t * }     % between first and second author (if more than two)
                  if$
                }
            if$
          }
          { t "\bibnf {others} {} {} {} {}" =
%            t "\bibnf {others} {} {} {} {} \editorsname" = or
            numnames #3 > or
               { "\bibimfont {\Bibetal}" * }
               { numnames #3 <
                   { "\Bibbtasep " * t * }  % between (exactly) two authors
                   { "\Bibbstasep " * t * } % between second and third author
                 if$
               }
            if$
          }
        if$
      }
     't
    if$
    nameptr #1 + 'nameptr :=
    namesleft #1 -  'namesleft :=
  }
  while$
}

FUNCTION {format.editor.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
  { s nameptr "{ll}" format.name$ 't := t editorlastnamemarkup 't :=
    s nameptr "{ff}" format.name$ othermarkup 'u :=
    s nameptr "{f.}" format.name$ othermarkup 'v :=
    s nameptr "{vv}" format.name$ othermarkup 'w :=
    s nameptr "{jj}" format.name$ othermarkup 'x :=
    t u * v * w * x * 't :=
     nameptr #1 >
      { namesleft #1 >
          { numnames #3 >
                { "" * }
                { numnames #3 <
                     { "\Bibbtasep " * t * }      % between (exactly) two authors
                     { "\Bibbfsasep " * t * }     % between first and second author (if more than two)
                  if$
                }
            if$
          }
          { t "\bibenf {others} {} {} {} {}" =
            t "\biblenf {others} {} {} {} {}" = or
            t "\bibenf {others} {} {} {} {} \editorsname" = or
            t "\biblenf {others} {} {} {} {} \editorsname" = or
            numnames #3 > or
               { "\bibimfont {\Bibetal}" * }
               { numnames #3 <
                   { "\Bibbtasep " * t * }        % between (exactly) two authors
                   { "\Bibbstasep " * t * }       % between second and third author
                 if$
               }
            if$
          }
        if$
      }
     't
    if$
    nameptr #1 + 'nameptr :=
    namesleft #1 -  'namesleft :=
  }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { pseudonym empty$                                                           % added JB
        { "" }                                                                   % added JB NOT 'skip$ !!!!
%        { "[" pseudonym * "]" * }                                               % added JB
         { "[" pseudonym format.names * "]" * }                                  % added JB  ?????
      if$                                                                        % added JB
    }
    { author format.names                                                        % added JB
      pseudonym empty$                                                           % added JB
        'skip$                                                                   % added JB NOT { "" } !!!!
%        { " [" * pseudonym * "]" * }                                            % added JB
        { " [" * pseudonym format.names * "]" * }                                % added JB  ????
      if$      
    }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { pseudonym empty$                                                           % added JB
         { "" }                                                                  % added JB NOT 'skip$ !!!!
%        { "[" pseudonym * "]" * }                                               % added JB
         { "[" pseudonym format.editor.names * "]" * }                           % added JB  ?????        
      if$                                                                        % added JB
    }
    { editor format.editor.names
      pseudonym empty$                                                           % added JB
        'skip$                                                                   % added JB NOT { "" } !!!!
%        { " [" * pseudonym * "]" * }                                            % added JB
         { " [" * pseudonym format.editor.names * "]" * }                        % added JB  ?????        
      if$                                                                        % added JB
      noed missing$       
          { editor num.names$ #1 >
                 { editortype empty$
                      { "\bibedformat {\editorsname}" * }
                      { "\bibedformat {" * editortype * "}" * } 
                   if$
                 }
                 { editortype empty$
                      { "\bibedformat {\editorname}" * }
                      { "\bibedformat {" * editortype * "}" * } 
                   if$
                 }
            if$
          }
          { noed 's :=
            s "1" =
                 'skip$
                 { editor num.names$ #1 >
                        { editortype empty$
                             { "\bibedformat {\editorsname}" * }
                             { "\bibedformat {" * editortype * "}" * } 
                          if$
                        }
                        { editortype empty$
                             { "\bibedformat {\editorname}" * }
                             { "\bibedformat {" * editortype * "}" * } 
                          if$
                        }
                   if$                 
                 }
            if$
          }
      if$
    }
  if$
}

FUNCTION {format.in.editors}
{ editor empty$
    { pseudonym empty$                                                           % added JB
        'skip$                                                                   % added JB
%         { "[" pseudonym * "]" * }                                              % added JB
          { "[" pseudonym format.editor.names * "]" * }                          % added JB  ?????        
      if$                                                                        % added JB
    }
    { editor format.editor.names
      pseudonym empty$                                                           % added JB
        'skip$                                                                   % added JB
%        { " [" * pseudonym * "]" * }                                            % added JB
         { " [" * pseudonym format.editor.names * "]" * }                        % added JB  ?????        
      if$                                                                        % added JB
      noed missing$       
          { editor num.names$ #1 >
                 { editortype empty$
                      { " \bibedinformat {\editorsname}" * }
                      { " \bibedinformat {" * editortype * "}" * } 
                   if$
                 }
                 { editortype empty$
                      { " \bibedinformat {\editorname}" * }
                      { " \bibedinformat {" * editortype * "}" * } 
                   if$
                 }
            if$
          }
          { noed 's :=
            s "1" =
                 'skip$
                 { editor num.names$ #1 >
                        { editortype empty$
                             { " \bibedinformat {\editorsname}" * }
                             { " \bibedinformat {" * editortype * "}" * } 
                          if$
                        }
                        { editortype empty$
                             { " \bibedinformat {\editorname}" * }
                             { " \bibedinformat {" * editortype * "}" * } 
                          if$
                        }
                   if$                 
                 }
            if$
          }
      if$
    }
  if$
}

% added JB
FUNCTION {format.titleaddon}
{ title empty$
    'skip$
    { titleaddon empty$
         'skip$
         { ", " * titleaddon * }
     if$
    }
  if$  
}

% added JB
FUNCTION {format.booktitleaddon}
{ booktitle empty$
    'skip$
    { booktitleaddon empty$
        'skip$
         { ", " * booktitleaddon * }
     if$
    }
  if$
}

FUNCTION {format.title}
{ before.all 'output.state := % added JB
  title empty$
    { "" }
    { title checked.emph.format swap$ * non.stop 
% JB, ONLY FOR CMS 
%       { before.all 'output.state := "\bibatsep} " *  }
%       { before.all 'output.state := "\bibatsep} " * }       
% CH 0.1g, superfluous \bibatsep  removed, reported by Maarten Wisse
       { before.all 'output.state := "}" * }
       { before.all 'output.state := "}" * }       
%        { "}" * }        
%        { "}" * }
      if$
    }
  if$
}

% added JB
FUNCTION {journalnamemarkup}
{ duplicate$ empty$
    { pop$ "" }
    { "\bibJTsep \bibjtfont {" swap$ * "}" * }
  if$
}

FUNCTION {seriesformat}
{ duplicate$ empty$
    { pop$ "" }
    { "\bibsnfont {" swap$ * "}" * }
  if$
}

FUNCTION {is.num}
{ chr.to.int$
  duplicate$ "0" chr.to.int$ < not
  swap$ "9" chr.to.int$ > not and
}

FUNCTION {extract.num}
{ duplicate$ 't :=
  "" 's :=
  { t empty$ not }
  { t #1 #1 substring$
    t #2 global.max$ substring$ 't :=
    duplicate$ is.num
      { s swap$ * 's := }
      { pop$ "" 't := }
    if$
  }
  while$
  s empty$
    'skip$
    { pop$ s }
  if$
}

FUNCTION {convert.edition.superscript}
{ edition extract.num "l" change.case$ 's :=
  s "first" = s "1" = or s "erste" = or
    { "1" 't := }
    { s "second" = s "2" = or s "zweite" = or
        { "2" 't := }
        { s "third" = s "3" = or s "dritte" = or
            { "3" 't := }
            { s "fourth" = s "4" = or s "vierte" = or
                { "4" 't := }
                { s "fifth" = s "5" = or s "fuenfte" = or
                    { "5" 't := }
                    { s "sixth" = s "6" = or s "sechste" = or
                        { "6" 't := }
                        { s "seventh" = s "7" = or s "siebte" = or s "siebente" = or
                            { "7" 't := }
                            { s "eighth" = s "8" = or s "achte" = or
                                { "8" 't := }
                                { s "nineth" = s "9" = or s "neunte" = or
                                      { "9" 't := }
                                      {  s #1 #1 substring$ is.num
                                           { s "" * 't := }
                                           { edition 't := }
                                         if$
                                      }
                                  if$
                                }
                              if$
                            }
                          if$
                        }
                      if$
                    }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
  t
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {output.bibitem.fst}
{ newline$
  "\bibitem[{" write$
  shortauthor empty$
    { juraauthor empty$
       { author empty$
          { editor empty$
              { organization empty$
                   { "empty" write$ }
                   { organization write$ }
                if$
              }
              { label write$ }
            if$
          }
          { label write$ }
         if$
       }
       { juraauthor write$ }
     if$
    }
    { shortauthor write$ }
  if$
  "}%" write$
  sortkey empty$
    'skip$
     { newline$ "% SORTKEY: [" sortkey * "]" * write$ }
  if$
  newline$
  shorttitle empty$
    { juratitle empty$
        { journal empty$
            { "         {}%" write$ }
            { title empty$
                  { key empty$
                       { "         {}%" write$ }
                       { "         {" key * "}%" * write$ }
                    if$
                  }
                  { "         {" title * "}%" * write$ }
              if$
            }
          if$
        }
        { "         {" juratitle * "}%" * write$ }
      if$
    }
    { "         {" shorttitle * "}%" * write$ }
  if$
  newline$
  repeat-authors
    { "         {{1}" write$ }
    { "         {{0}" write$ }
  if$
  howcited empty$
    { "{}" write$ }
    { "{" howcited * "}" * write$ }
  if$
  type$ empty$
    { "{}" write$ }
    { "{" type$ * "}" * write$ }
  if$
  year empty$
    { "{}" write$ }
    { "{" year * write$
      extra.label "" =
          'skip$
          { "\el {" write$ extra.label write$ "}" write$ }
      if$
      "}" write$
    }
  if$  
  edition empty$
    { "{}" write$ }
    { "{" convert.edition.superscript extract.num * "}" * write$ }
  if$
  ssedition empty$
    { "{}" write$ }
    { "{" ssedition * "}" * write$ }
  if$
  language empty$
    { "{}" write$ }
    { "{" language * "}" * write$ }
  if$
  pages empty$
    { "{}%" write$ }
    { "{" write$ pages n.dashify write$ "}%" write$ }
  if$
  newline$
}

FUNCTION {output.bibitem.snd}
{ publisher empty$
    { address empty$
         { year empty$
               { "          {" write$ }
               { "          {" write$ 
                            year write$ 
                            extra.label "" =
                              'skip$
                              { "\el {" write$ extra.label write$ "}" write$ }
                            if$                               
               }
           if$
         }
         { year empty$
               { "          {" write$ address write$ }
               { "          {" write$ address write$ "\bibbdsep {} " write$
                            year write$
                            extra.label "" =
                              'skip$
                              { "\el {" write$ extra.label write$ "}" write$ }
                            if$               
               }
           if$
         }
      if$
    }
    { address empty$
        { year empty$
              { "          {" write$ publisher write$ }
              { "          {" write$ publisher write$ "\bibbdsep {} " write$ 
                            year write$
                            extra.label "" =
                              'skip$
                              { "\el {" write$ extra.label write$ "}" write$ }
                            if$
              }
          if$
        }
        { year empty$
              { "          {" write$ address write$ "\bpubaddr {} " write$ publisher write$ }
              { "          {" write$ address write$ "\bpubaddr {} " write$ publisher write$ "\bibbdsep {} " write$ 
                            year write$
                            extra.label "" =
                              'skip$
                              { "\el {" write$ extra.label write$ "}" write$ }
                            if$
              }
          if$
        }
      if$
    }
  if$
  "}}%" write$ newline$
  "         {" write$
  title empty$
    { "{}%" write$ newline$ }
    { "{" title * "}%" * write$ newline$ }
  if$
  gender empty$
    { "          {}" write$ }
    { "          {" gender * "}" * write$ }
  if$
  "{}{}{}{}{}{}{}}%" write$ newline$
  "        ]{" write$
  cite$ write$
  "}" write$
  url empty$
    'skip$
    { newline$ " \jburldef {" write$ cite$ write$ "}{" write$ url write$ "}%" write$ }
  if$
  newline$
  " \jbbibargs {" write$
  ""
  before.all 'output.state :=
}

FUNCTION {output.bibitem}{ output.bibitem.fst output.bibitem.snd }

FUNCTION {word.in}
{ bbl.in
  " " * }

% JB extra label stuff added in 0.1a
FUNCTION {add.extra.label}
{extra.label "" =
   'skip$
   { "\bibel {" extra.label * "}" * * }
 if$
}

FUNCTION {format.date}
{ year empty$
    'skip$
    { output.state before.all =                                                          % added JB to avoid leading whitespace for minimal misc
            { "" * mid.sentence 'output.state := }                                       % added JB to avoid leading whitespace for minimal misc
            { " " * }                                                                    % added JB to avoid leading whitespace for minimal misc
      if$                                                                                % added JB to avoid leading whitespace for minimal misc
      "\dateldelim {}" *
      month empty$
          'skip$
          { month "{} " * * }
      if$
      year * add.extra.label "\daterdelim {}" * 
    }
  if$
}

FUNCTION {format.btitle}
{ before.all 'output.state :=                                                            % added JB to avoid leading comma
  title emphasize
}

% NEW 0.0.1h 
FUNCTION {format.booktitle}
{ before.all 'output.state :=                                                            % added JB to avoid leading comma
  booktitle b.emphasize
}


FUNCTION {format.booklet.title}
{ before.all 'output.state :=                                                            % added JB to avoid leading comma
  title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

% added JB
FUNCTION {format.bvolume}
{ volume empty$
    { " " }
    { 
      volumetitle empty$
        { 
%           type$ "incollection" = type$ "inbook" = or % Added JB 0.0.1f, not sure about the types ...
%              { volume }                              % no `volume' for missing volumetitle
%              {                                       % Removed again in 0.0.1m (???)
%                NEW 0.0.1n no \volumename for crossreff'ing works
                 crossref missing$
                     { bbl.volume volume tie.or.space.connect }
% NEW 0.0.1o no bbl.vol=volname added
                     { bbl.vol volume tie.or.space.connect }
                 if$
%              } 
%          if$
        }
        {
%                NEW 0.0.1n no \volumename for crossreff'ing works         
          crossref missing$
              { bbl.volume volume tie.or.space.connect }
% NEW 0.0.1o no bbl.vol=volname added
              { bbl.vol volume tie.or.space.connect }
          if$
          ": \bibvtfont {" * volumetitle "}" * * 
        }
      if$
% JB removed       "volume and number" number either.or.check
    }
  if$
}

% added JB
FUNCTION {format.bvolume.btitle.reversed}
{ volume empty$
    { title empty$
        'skip$
        { 
          before.all 'output.state :=
          title emphasize 
        }
      if$    
    }
    { volumetitle empty$
        { title empty$
             'skip$
             { 
                before.all 'output.state := 
                title emphasize
             }
          if$
        }
        { "\bibvtfont {" volumetitle * "}" * *
          bbl.volume volume tie.or.space.connect 
          title empty$
            'skip$
            { bbl.of title emphasize * * }
          if$
        }
      if$
    }
  if$
}

% added for volformat handling by JB, for @book and @inbook
FUNCTION {format.bvolume.btitle}
{ volformat missing$
     { 
       format.btitle "title" output.check format.bvolume output
     }
     { volformat 's := 
       s "reversed" =
            { volumetitle empty$
                 { 
                   "volformat given but no volumetitle in " cite$ * warning$ 
                   format.bvolume.btitle.reversed output
                 }
                 { 
                   volume empty$
                        { "there's a volumetitle but no volume in " cite$ * warning$ * }
                        'skip$                 
                   if$
                   title empty$
                     { 
                        "volformat given but no title in " cite$ * warning$ 
                     } 
                     { 
                        format.bvolume.btitle.reversed output
                     }
                   if$
                 }
              if$
            }
            { 
                format.btitle "title" output.check format.bvolume output
            }            
       if$
     }     
  if$
}


% added JB, for @incollection (booktitle!)
FUNCTION {format.bvolume.booktitle.reversed}
{ volume empty$
    { booktitle empty$
        'skip$
        { 
          before.all 'output.state :=
          booktitle b.emphasize 
        }
      if$    
    }
    { volumetitle empty$
        { booktitle empty$
             'skip$
             { 
                before.all 'output.state := 
                bbl.volume volume tie.or.space.connect
                booktitle empty$
                    'skip$
                    { bbl.of booktitle b.emphasize * * }
                if$                
             }
          if$
        }
        { "\bibvtfont {" volumetitle * "}" * *
          bbl.volume volume tie.or.space.connect 
          booktitle empty$
            'skip$
            { bbl.of booktitle b.emphasize * * }
          if$
        }
      if$
    }
  if$
}

% added for volformat handling by JB
FUNCTION {format.in.editor.booktitle.volume.voltitle}
{ editor empty$
      { chapter empty$                                                              % no comma between chap. and `in:'
           { ", " * }
           { " " *  }
        if$
        word.in *
      }
      { chapter empty$
           { ", " * }
           { " " *  }
        if$
        word.in format.in.editors * ", " * * 
      }
  if$
  volformat missing$
     { 
       format.booktitle "booktitle" output.check
       format.bvolume output                                                         % 
       format.booktitleaddon                                                         % added JB 0.0.1m booktitleaddon
     }                                                                               % after volume
     { volformat 's := 
       s "reversed" =
            { volumetitle empty$
                 { 
                    "volformat given but no volumetitle in " cite$ * warning$ 
                    format.bvolume.booktitle.reversed output
                    format.booktitleaddon                                            % added JB 0.0.1m                      
                 }
                 { 
                   volume empty$
                        { "there's a volumetitle but no volume in " cite$ * warning$ * }
                        'skip$                 
                   if$
                   booktitle empty$
                     { 
                        "volformat given but no booktitle in " cite$ * warning$ 
                     } 
                     { 
                        format.bvolume.booktitle.reversed output 
                        format.booktitleaddon                                        % added JB 0.0.1m                         
                     }
                   if$
                 }
              if$
            }
            { 
               format.booktitle "booktitle" output.check
               format.bvolume output                                                 %
               format.booktitleaddon                                                 % added JB 0.0.1m booktitleaddon
            }                                                                        % after volume
       if$
     }     
 if$
}

% new function by JB
FUNCTION {format.number.series.book}
{ number empty$
   { series seriesformat field.or.null }        
   { series empty$
       { "there's a number but no series in " cite$ * warning$ number }
       { series seriesformat number tie.or.space.connect }
     if$
   }
 if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
% JB removed         { series field.or.null }
        { series seriesformat field.or.null }        
%                ^^^^^^^^^^^^ added JB
        { output.state mid.sentence =
            { bbl.number }
            { bbl.number capitalize }
          if$
          number tie.or.space.connect
          series empty$
            { "there's a number but no series in " cite$ * warning$ }
% JB removed             { bbl.in space.word * series * }
            { bbl.in.series space.word * series seriesformat * }
%                                               ^^^^^^^^^^^^ added JB            
          if$
        }
      if$
    }
    { type$ "incollection" = crossref missing$ and
        { format.number.series.book }                                               % to make output equal to crossref ...
        { "" }
      if$ 
   }
   if$
}

FUNCTION {first.ed} { "\firstedname" }
FUNCTION {second.ed}{ "\secondedname" }
FUNCTION {third.ed} { "\thirdedname" }
FUNCTION {fourth.ed}{ "\fourthedname" }
FUNCTION {fifth.ed} { "\fifthedname" }
FUNCTION {th.ed}    { "\thedname" }

FUNCTION {convert.edition}
{ edition extract.num "l" change.case$ 's :=
  s "first" = s "1" = or
  { first.ed 't := }
  { s "second" = s "2" = or
    { second.ed 't := }
    { s "third" = s "3" = or
      { third.ed 't := }
      { s "fourth" = s "4" = or
        { fourth.ed 't := }
        { s "fifth" = s "5" = or
          { fifth.ed 't := }
          { s #1 #1 substring$ is.num
            { s th.ed * 't := }
            { edition 't := }
            if$
          }
          if$
        }
        if$
      }
      if$
    }
    if$
  }
  if$
  t
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { convert.edition "~\editionname {}" * }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { before.all 'output.state := % added JB
% JB NO space before \jbPage(s), beause of a additional space after 
% \@gobbl'ing their in citations, space is added by definition of \jbPage(s)!!!
      "\jbPages {" pages n.dashify * "}" * 
    }
  if$
}

FUNCTION {format.journal.pages}
{ pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { " \jbartPages {" pages n.dashify * "}" * * }
      if$
    }
  if$
}

% added JB, pages after year!
% \artvolnumformat introduced
FUNCTION {format.vol.num.pages}
{ before.all 'output.state :=
  volume empty$
    { year empty$ not
         { "there's a year but no volume in " cite$ * warning$ }
         'skip$
      if$
      number empty$
          { "" }
          { "there's a number but no volume in " cite$ * warning$ 
            "~\artvolnumformat {}{" number * "}" * }
      if$
    }
    { number empty$
          { "~\artvolnumformat {" volume * "}{}" * }
          { "~\artvolnumformat {" volume * "}{" number * "}" * * }
      if$
    }
  if$
  year empty$
    { month empty$ not
         { "there's a month but no year in " cite$ * warning$ }
         'skip$
      if$
    }
    { month empty$
        { " \artyearformat {" year * "}" * * }
        { " \artyearformat {" month * "{} " * year * "}" * * }
      if$
    }
  if$
% JB removed to abe able to put url before pages
%  format.journal.pages
}

% new JB special function to place only chaptern not chapter and pages before crossref stuff for inbook
FUNCTION {format.chapter.inbook.incoll}
{ chapter empty$
    { "" }
    { before.all 'output.state :=
      type empty$
        { ", " bbl.chapter * }
        { ", " type * "l" change.case$ }
      if$
      chapter tie.or.space.connect
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
% JB removed (INBOOK)       { bbl.chapter }    
        { before.all 'output.state := ", " bbl.chapter * }
% JB removed (INBOOK)       { type "l" change.case$ }
        { ", " type * "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
        'skip$
% JB removed         { ", " * format.pages * }        
        { format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
        { word.in booktitle b.emphasize format.booktitleaddon * }                            % booktitleaddon added in 0.0.1m
        { word.in format.in.editors * ", " * booktitle b.emphasize format.booktitleaddon * } % booktitleaddon added in 0.0.1m
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { bbl.techrep }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

% added JB
FUNCTION {format.article.crossref}
{ " \protect\jbartcrossrefchecked {" crossref * "}" * }

% changed JB
FUNCTION {format.book.crossref}% only used by @book and @inbook
{ type$ "book" = % make this difference (\jbcrossrefchecked[<volume>]{<key>}) only for books
    { volume empty$
        { "empty volume in " cite$ * "'s crossref of " * crossref * warning$ 
% JB removed     word.in
          "\protect\jbcrossrefchecked {" crossref * "}" * }
        { "\protect\jbcrossrefchecked[" volume * "]{" * crossref "}" * * }
      if$
    }
    { "\protect\jbcrossrefchecked {" crossref * "}" * }
  if$
}

% added JB
%FUNCTION {format.incoll.inproc.crossref.volume}
%{ " \protect\jbincollcrossref {" crossref * "}" * 
%  % check added, if editor for the crossref'd work is given, switch introduced for that problem.
%  editor empty$
%    { "{}" * }
%    { "{1}" * }
%  if$
% "{" * volume * "}{" * volumetitle * "}" *
%}

FUNCTION {format.incoll.inproc.crossref}
{ " \protect\jbincollcrossref {" crossref * "}" * 
  % check added, if editor for the crossref'd work is given, switch introduced for that problem.
  editor empty$
    { "{}" * }
    { "{1}" * }
  if$
}

FUNCTION {oaddress.opublisher.oyear.format}
{ type$ "book" = 
% ???? type$ "inbook" = type$ "incollection" = or or
     { oyear empty$
            { opublisher empty$
                 { oaddress empty$
                     'skip$
                     { oaddress "\osep {} " * * }
                   if$
                 }
                 { oaddress empty$
                     { opublisher "\osep {} " * * }
                     { oaddress "\bpubaddr {} " * opublisher "\osep {} " * * * }
                   if$
                 }
              if$
            }
            { year empty$
                 { "oyear without year in " cite$ * warning$ }
                 'skip$
              if$
              opublisher empty$
                 { oaddress empty$
                     { oyear * "\osep {} \reprintname {}, " * }
                     { oaddress "\bibbdsep {} " * oyear "\osep {} " * * * }
                   if$
                 }
                 { oaddress empty$
                     { opublisher * "\bibbdsep {} " * oyear "\osep {} " * * }
                     { oaddress "\bpubaddr {} " * opublisher "\bibbdsep {} " * * oyear "\osep {} " * * * }
                   if$
                 }
              if$
            }
       if$
     }
     'skip$
  if$
}

% added JB (completely rewritten function, cares about month; somewhat simple, but better to read)
% \el stuff added in v0.1a
FUNCTION {format.publisher}
{ year empty$ month empty$ not and 
        { "there's a month but no year in " cite$ * warning$ }
        'skip$
  if$
  address empty$ publisher empty$ year empty$ and and
     'skip$
     { " \apyformat {" * }
  if$
  oaddress.opublisher.oyear.format                                              
  address empty$
    { 
      publisher empty$
          {
            year empty$
                'skip$
                { month empty$
                     { year * add.extra.label }
                     { month "{} " * year * * add.extra.label }
                  if$
                }
            if$
          }
          {
            year empty$
                { month empty$
                     { publisher * }
                      'skip$                                                  
                  if$
                }
                { month empty$
                     { publisher "\bibbdsep {} " * year * * add.extra.label }
                     { publisher "\bibbdsep {} " * month "{} " * year * * * add.extra.label }
                  if$
                }                    
            if$
          }
      if$
    }
    { 
      publisher empty$
          { 
            year empty$
                { address * }
                { month empty$
                     { address "\bibbdsep {} " * year * * add.extra.label }
                     { address "\bibbdsep {} " * month "{} " * year * * * add.extra.label }
                  if$
                }
            if$
          }
          { 
            year empty$
                { month empty$
                     { address "\bpubaddr {} " * publisher * * }
                      'skip$                                                  
                  if$
                }
                { month empty$
                     { address "\bpubaddr {} " * publisher "\bibbdsep {} " * year * * * add.extra.label }
                     { address "\bpubaddr {} " * publisher "\bibbdsep {} " * month "{} " * year * * * * add.extra.label }
                  if$
                }                    
            if$
          }
      if$
    }        
  if$
  address empty$ publisher empty$ year empty$ and and
     'skip$
     { "}" * }
  if$
}

% added JB (completely rewritten function, cares about month; somewhat simple, but better to read)
% TODO: try to make one function for all (..., manual, techreport)
FUNCTION {format.publisher.manual}
{ year empty$ month empty$ not and 
        { "there's a month but no year in " cite$ * warning$ }
        'skip$
  if$
  address empty$ organization empty$ year empty$ and and
     'skip$
     { " \apyformat {" * }
  if$
  address empty$
    { 
      organization empty$
          {
            year empty$
                'skip$
                { month empty$
                     { year * add.extra.label }
                     { month "{} " * year * * add.extra.label }
                  if$
                }
            if$
          }
          {
            year empty$
                { month empty$
                     { organization * }
                      'skip$                                                  
                  if$
                }
                { month empty$
                     { organization "\bibbdsep {} " * year * * add.extra.label }
                     { organization "\bibbdsep {} " * month "{} " * year * * * add.extra.label }
                  if$
                }                    
            if$
          }
      if$
    }
    { organization empty$
          { year empty$
                { address * }
                { month empty$
                     { address "\bibbdsep {} " * year * * add.extra.label }
                     { address "\bpubaddr {} " * month "{} " * year * * * add.extra.label }
                  if$
                }
            if$
          }
          { year empty$
                { month empty$
                     { address "\bpubaddr {} " * organization * * }
                      'skip$                                                  
                  if$
                }
                { month empty$
                     { address "\bpubaddr {} " * organization "\bibbdsep {} " * year * * * add.extra.label }
                     { address "\bpubaddr {} " * organization "\bibbdsep {} " * month "{} " * year * * * * add.extra.label }
                  if$
                }                    
            if$
          }
      if$
    }        
  if$
  address empty$ organization empty$ year empty$ and and
     'skip$
     { "}" * }
  if$
}

% added JB (completely rewritten function, cares about month; somewhat simple, but better to read)
% TODO: try to make one function for all (..., manual, techreport)
FUNCTION {format.publisher.techreport}
{ year empty$ month empty$ not and 
        { "there's a month but no year in " cite$ * warning$ }
        'skip$
  if$
  address empty$ institution empty$ year empty$ and and
     'skip$
     { " \apyformat {" * }
  if$
  address empty$
    { 
      institution empty$
          {
            year empty$
                'skip$
                { month empty$
                     { year * add.extra.label }
                     { month "{} " * year * * add.extra.label }
                  if$
                }
            if$
          }
          {
            year empty$
                { month empty$
                     { institution * }
                      'skip$                                                  
                  if$
                }
                { month empty$
                     { institution "\bibbdsep {} " * year * * add.extra.label }
                     { institution "\bibbdsep {} " * month "{} " * year * * * add.extra.label }
                  if$
                }                    
            if$
          }
      if$
    }
    { institution empty$
          { year empty$
                { address * }
                { month empty$
                     { address "\bibbdsep {} " * year * * add.extra.label }
                     { address "\bpubaddr {} " * month "{} " * year * * * add.extra.label }
                  if$
                }
            if$
          }
          { year empty$
                { month empty$
                     { address "\bpubaddr {} " * institution * * }
                      'skip$                                                  
                  if$
                }
                { month empty$
                     { address "\bpubaddr {} " * institution "\bibbdsep {} " * year * * * add.extra.label }
                     { address "\bpubaddr {} " * institution "\bibbdsep {} " * month "{} " * year * * * * add.extra.label }
                  if$
                }                    
            if$
          }
      if$
    }        
  if$
  address empty$ institution empty$ year empty$ and and
     'skip$
     { "}" * }
  if$
}


% added JB
FUNCTION {format.edby.editors}
{ author empty$                                                                     
    'skip$                                                                          
    { editor empty$                                                                 
         'skip$                                                                     
         { "\edbyname {} " editor format.editor.names * output }
      if$                                                                           
    }                                                                               
  if$                                                                               
}

FUNCTION {make.translators}
{ 
% JB removed   before.all 'output.state :=
  translator empty$
    { flanguage empty$
        { "" }
        { "\translator {}{" flanguage "l" change.case$ * "}{1}" * }
      if$
    }
    { 
%      "\translator {" translator format.names * "}" *
      "\translator {" translator * "}" *      
      flanguage empty$
        { "{}{1}" * }
        { "{" * flanguage "l" change.case$ * "}{1}" * }
      if$
    }
  if$
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name
  format.title "title" output.check
  format.titleaddon                                                                 % added JB because of moving format.titleaddon outside format.(b)title
  crossref missing$
    { journal journalnamemarkup                                                     % added JB
      "journal" output.check
      format.vol.num.pages output
      write.url                                                                     % added JB
      format.journal.pages                                                          % added JB (extracted from format.vol.num.pages to place url before pages)
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  format.note                                                                       % added JB
  fin.entry
}

FUNCTION {review}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name
  format.title "title" output.check
  format.titleaddon                                                                 % added JB because of moving format.titleaddon outside format.(b)title
  crossref missing$
    { journal journalnamemarkup                                                     % added JB
      "journal" output.check
      format.vol.num.pages output
      write.url                                                                     % added JB
      format.journal.pages                                                          % added JB (extracted from format.vol.num.pages to place url before pages)
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  format.note                                                                       % added JB
  fin.entry
}

FUNCTION {format.editors.last}
{ editor empty$
   { "" before.all 'output.state := }
%       ^^^^^^^^^^^^^^^^^^^^^^^^ CH 0.1j
   { author empty$ not
        { " \edbyname{} " editor format.editor.names * }
        { "" before.all 'output.state := }                                        % CH 0.1h from 'skip$ to { "" }!
%            ^^^^^^^^^^^^^^^^^^^^^^^^ CH 0.1j
    if$
   }
  if$
}

FUNCTION {book}
{ output.bibitem
  author empty$
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CH 0.1i %%%%%%%%%%%%%%%%%%%%%%%%                   % addedd JB, editors should appear first, if authors are missing 
    { editor empty$ not
	 { format.editors output.nonnull                                           % JB output.checked removed 
         }                                                                         
	 { "neither author nor editor in " cite$ * warning$ }                      % added JB, print out warning only if neither author nor editor given
      if$   
    }
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CH 0.1i %%%%%%%%%%%%%%%%%%%%%%%%
    { format.authors output.nonnull }
  if$
  set.after.author.sep output                                                       % added JB to set "}{" after the name
  crossref missing$
    { format.bvolume.btitle                                                         % added JB, complete handling of volume and volumetitle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CH 0.1i %%%%%%%%%%%%%%%%%%%%%%%%%%  
% author empty$ editor empty$ not and { 					    % removed JB, checks are done inside of format.editors.last 
    format.editors.last output.nonnull                                              % JB output.checked removed 
%    } 'skip$ if$                                                                   % removed JB, checks are done inside of format.editors.last
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CH 0.1i %%%%%%%%%%%%%%%%%%%%%%%%%
      make.translators output                                                       % added JB      
      format.titleaddon                                                             % added JB      
      format.number.series.book output                                              % added JB
      format.edition output
% added in 0.1b
      format.pages output
      format.publisher                                                              % added JB 
      write.url                                                                     % added JB (url before note) 
      format.note                                                                   % added JB  
    }
    {
      format.btitle "title" output.check
      format.book.crossref output.nonnull
    }
  if$
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.btitle "title" output.check
  format.titleaddon                                                                 % added JB because of moving format.titleaddon outside format.(b)title
  howpublished output
  format.publisher
% added in 0.1b
  format.pages output  
  write.url                                                                         % added JB (url before note)
  format.note                                                                       % added JB  
  fin.entry
}

% JB removed FUNCTION {inbook}
% JB removed { output.bibitem
% JB removed   author empty$
% JB removed     { format.editors "author and editor" output.check
% JB removed % JB removed       editor format.key output
% JB removed     }
% JB removed     { format.authors output.nonnull
% JB removed       crossref missing$
% JB removed         { "author and editor" editor either.or.check }
% JB removed         'skip$
% JB removed       if$
% JB removed     }
% JB removed   if$
% JB removed   set.after.author.sep output                                                       % added JB to set "}{" after the name  
% JB removed   format.btitle "title" output.check
% JB removed   crossref missing$
% JB removed     {
% JB removed       format.bvolume output
% JB removed % JB removed       format.chapter.pages "chapter and pages" output.check            % chapter and pages together before publisher????
% JB removed       format.chapter.inbook.incoll output                                           % put out only chapter before crossref stuff ...
% JB removed       format.number.series output
% JB removed       format.edition output
% JB removed       format.titleaddon                                                             % added JB
% JB removed % JB removed       note output
% JB removed       format.publisher
% JB removed       format.pages output                                                           % ... and put the pages after publisher      
% JB removed       write.url                                                                     % added JB, moved up into non-crossref way othwerwise double url's
% JB removed       format.note                                                                   % added JB        
% JB removed     }
% JB removed     {
% JB removed % JB removed       format.chapter.pages "chapter and pages" output.check            %
% JB removed       format.titleaddon                                                             % added JB
% JB removed       format.chapter.inbook.incoll output                                           % put out only chapter before crossref stuff ...
% JB removed       format.book.crossref output.nonnull
% JB removed       format.pages output                                                           % ... and put the pages after crossref 
% JB removed     }
% JB removed   if$
% JB removed %  note output %note moved
% JB removed % JB removed   write.url
% JB removed   fin.entry
% JB removed % JB removed   write.url  
% JB removed }

FUNCTION {inbook} % copy of incollection v0.0.1d, only diff is the formatting of the title 
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  crossref missing$
    { format.bvolume.btitle                                                         % added JB, complete handling of volume and volumetitle
%      format.in.ed.booktitle "booktitle" output.check                              % removed, no need for booktitle
%      format.bvolume output
%      format.booktitleaddon
      format.titleaddon                                                             % added JB
      make.translators output                                                       % added JB
      format.number.series output
      format.edition output
      format.publisher
      write.url                                                                     % added JB, moved up into non-crossref way othwerwise double url's
      format.chapter.inbook.incoll output                                           % chapter moved back to pages
      format.pages output                                                           %     
      format.note                                                                   % added JB      
    }
    { format.btitle "title" output.check
%            ^ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      format.chapter.inbook.incoll output                                           % put out only chapter before crossref stuff ...
      format.incoll.inproc.crossref output.nonnull
      format.bvolume output                                                         % <volume> and <voluemtitle> after crossref      
%      format.incoll.inproc.crossref.volume output.nonnull                           % <volume> and <voluemtitle> inside crossref ?
      format.pages output                                                           % ... and put the pages after crossref       
    }
  if$
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  crossref missing$
    { format.btitle "title" output.check
      format.titleaddon                                                             % ????? because of removing format.titleaddon from format.title
                                                                                    % and placement inside crossref missing$-way.
      format.chapter.inbook.incoll output                                           % put out only chapter before `in:' ...
      format.in.editor.booktitle.volume.voltitle                                    % added JB, complete handling of volume and volumetitle
      make.translators output                                                       % added JB      
      format.number.series.book output
      format.edition output
% JB removed       format.chapter.pages output
      format.publisher
      write.url                                                                     % added JB, moved up into non-crossref way othwerwise double url's
      format.note                                                                   % added JB      
      format.pages output                                                           % ... and put the pages after all      
    }
    { format.title "title" output.check
      format.chapter.inbook.incoll output                                           % put out only chapter before crossref stuff ...
      format.incoll.inproc.crossref output.nonnull
      format.bvolume output                                                         % <volume> and <voluemtitle> after crossref      
%      format.incoll.inproc.crossref.volume output.nonnull                           % <volume> and <voluemtitle> inside crossref ?
      format.pages output                                                           % ... and put the pages after crossref       
% JB removed       format.chapter.pages output
    }
  if$
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.title "title" output.check
  format.titleaddon                                                                 % ????? because of removing format.titleaddon from format.title  
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
% JB      format.pages output
% JB      note output
      publisher empty$
        { organization output
          format.publisher
        }
        { organization output
          format.publisher
        }
      if$
      write.url                                                                     % added JB
      format.note                                                                   % added JB
      format.pages output                                                           % added JB      
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
%  note output % note moved
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  format.authors output
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.btitle "title" output.check
  format.titleaddon                                                                 % added JB
  format.edition output
  format.publisher.manual                                                           % added JB
  write.url                                                                         % added JB
  format.note                                                                       % added JB
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.btitle "title" output.check
  format.titleaddon                                                                 % added JB
  bbl.mthesis format.thesis.type output.nonnull
  school "school" output.check
  format.publisher
  write.url                                                                         % added JB
  format.note                                                                       % added JB  
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.btitle output % JB removed "title" output.check                            % added JB
  format.titleaddon                                                                 % added JB
  howpublished output
  format.date % JB removed "year" output.check.special
  write.url                                                                         % added JB
  format.note                                                                       % added JB  
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.btitle "title" output.check
  format.titleaddon                                                                 % added JB  
  bbl.phdthesis format.thesis.type output.nonnull
  school "school" output.check
  format.publisher
  write.url                                                                         % added JB
  format.note                                                                       % added JB
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  format.editors output
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.btitle "title" output.check
  format.titleaddon                                                                 % added JB
% JB removed   format.bvolume output
% JB removed   format.number.series output
   format.number.series.book output                                                 % added JB to use the <series> <number> format
  organization output
  format.publisher
  write.url                                                                         % added JB
  format.note                                                                       % added JB
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.title "title" output.check
  format.tr.number output.nonnull
% JB removed   institution "institution" output.check
% JB removed   format.publisher
  format.publisher.techreport                                                       % added JB
  write.url                                                                         % added JB
  format.note                                                                       % added JB  
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  set.after.author.sep output                                                       % added JB to set "}{" after the name  
  format.title "title" output.check
  format.titleaddon                                                                 % added JB  
  format.date % JB removed "year" output.check.special
  write.url                                                                         % added JB
  format.note                                                                       % added JB
  fin.entry
}

% added JB
FUNCTION {www}
{ output.bibitem
  format.authors output.nonnull % "author" output.check % removed JB (to avoid warning about missing author)
  set.after.author.sep output
%  before.all 'output.state :=
  format.title output.nonnull   % "title" output.check  % removed JB (to avoid warning about missing title)
  write.url                                                                         % added JB  
  format.note                                                                       % added JB
  fin.entry
}

FUNCTION {default.type} { misc }

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}


% JB added, contributed by Oren Patashnik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INTEGERS { nextchar-ptr }                % points just past 1st "text" char

FUNCTION {get.nextchar.ptr}              % set nextchar-ptr to be just past
{ #1 'nextchar-ptr :=                    % the first "text" char---this will
    { duplicate$                         % be not quite trivial if
      nextchar-ptr #1 substring$         % the first "text" char is a
      "" = not                           % special characters like {\'E}
    }
    { nextchar-ptr #1 + 'nextchar-ptr := }
  while$
  pop$
  nextchar-ptr
}

STRINGS { first-char }                   % the to-be-capitalized character

FUNCTION {capitalize.von.part}           % capitalize 1st character of string
{ duplicate$ #1 text.prefix$             % get 1st character
  duplicate$ get.nextchar.ptr
  swap$ "u" change.case$                 % capitalize it
  'first-char :=
  global.max$ substring$                 % get rest of string
  first-char swap$ *                     % append the two
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% JB added
FUNCTION {format.lab.names}
{ 's :=                                     % for two: (Jones/de~Bruijn)
  s #1 "{vv~}{ll}" format.name$ capitalize.von.part % for three: (Jones/de~Bruijn/Smith)
  s num.names$ duplicate$                  % for more: (Jones et~al.)
  #2 <
  'pop$                                    % for just one person
  { duplicate$ duplicate$
    #2 >
    swap$ s swap$ "{ff }{vv }{ll}{ jj}" format.name$ "others" =
    or
    { pop$ "\Bibetal" * }                  % for 4 or more, or last is `others'
    { #2 =
      {  "\jbbtasep " * s #2 "{vv~}{ll}" format.name$ capitalize.von.part * } % for (exactly) two
      { "\jbbfsasep " * s #2 "{vv~}{ll}" format.name$ capitalize.von.part *   % for three
        "\jbbstasep " * s #3 "{vv~}{ll}" format.name$ capitalize.von.part *
      }
      if$
    }
    if$
  }
  if$
% NEW 0.1b  
  year missing$
      'skip$
     { "\jbdy {" year * "}" * * }% dummy year to prevent from `1993a' and `2001b' for works by the same author
  if$      
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            'key
          if$
        }
        { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.label}
{ type$ "book" = type$ "inbook" = or type$ "incollection" = or
  'author.editor.key.label
  { type$ "proceedings" =
    'editor.key.label
    'author.key.label
    if$
  }
  if$
  'label :=
}

% JB added in 0.1a
FUNCTION {calc.dummy.label}
{ type$ "book" = type$ "inbook" = or type$ "incollection" = or type$ "commented" = or
  'author.editor.key.label
  { type$ "proceedings" =
    'editor.key.label
    'author.key.label
    if$
  }
  if$
  'dummy.label :=
}

% JB added
FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  { namesleft #0 > }
  { nameptr #1 >
  { "   " * }
  'skip$
    if$
    s nameptr "{ll{ }}{  ff{ }}{ vv{ }}{  jj{ }}" format.name$ 't :=
    nameptr numnames = t "\bibnf {others} {} {} {} {}" = and
  { "\Bibetal" * }
  { t sortify * }
    if$
    nameptr #1 + 'nameptr :=
    namesleft #1 - 'namesleft :=
  }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

% added JB, pseudo code added
FUNCTION {author.sort}
{ author empty$
    { pseudonym empty$
          { key empty$
               { "to sort, need author, pseudonym or key in " cite$ * warning$ "" }
               { key sortify }
            if$
          }
          { pseudonym sort.format.names }
       if$
    }
    { author sort.format.names }
  if$
}

% added JB, pseudo code added
FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { pseudonym empty$
             { key empty$
                  { "to sort, need author, editor, pseudonym or key in " cite$ * warning$ "" }
                  { key sortify }
               if$
             }
             { pseudonym sort.format.names }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

% added JB, pseudo code added
FUNCTION {editor.sort}
{ editor empty$
    { pseudonym empty$
         { key empty$
              { "to sort, need author, editor, pseudonym or key in " cite$ * warning$ "" }
              { key sortify }
           if$
         }
         { pseudonym sort.format.names }
      if$
    }
    { editor sort.format.names }
  if$
}

% JB added
FUNCTION {presort}
{ calc.label
% JB added in 0.1a
  calc.dummy.label
  sortkey empty$
     { type$ "book" =
         type$ "incollection" =    or
           type$ "inbook" =        or
            'author.editor.sort
            { type$ "proceedings" =
                'editor.sort
                'author.sort
              if$
            }
        if$
     }
     { sortkey sort.format.names }
  if$
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

%FUNCTION {before.sort} { "BEFORE SORT:" top$ }
%EXECUTE {before.sort}
%FUNCTION {print.sort.keys} {sort.key$ top$}
%ITERATE {print.sort.keys}

SORT

%FUNCTION {after.sort} {"AFTER SORT:" top$}
%EXECUTE {after.sort}
%ITERATE {print.sort.keys}

%%|< JB added 
% JB added in 0.1a
STRINGS { last.dummy.label next.extra }
INTEGERS { last.extra.num number.label }
FUNCTION {forward.pass}
{ last.dummy.label dummy.label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      dummy.label 'last.dummy.label :=
    }
  if$
  number.label #1 + 'number.label :=
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  extra.label 'next.extra :=
  extra.label
  duplicate$ empty$
    'skip$
    { "" swap$ * }
  if$
  'extra.label :=
  dummy.label extra.label * 'dummy.label :=
}

ITERATE {forward.pass}
REVERSE {reverse.pass}


FUNCTION {set.current.authors}
{ author empty$
  { editor empty$
    { organization empty$
      { cite$ 'current-authors := }
      { organization 'current-authors := }
      if$
    }
    { editor format.names #1 entry.max$ substring$ 'current-authors := }
    if$
  }
  { author format.names #1 entry.max$ substring$ 'current-authors := }
  if$
}

STRINGS { previous-authors }

FUNCTION {initialize.previous.authors} { "" 'previous-authors := }

FUNCTION {mark.repeat.authors}
{ current-authors previous-authors =
  { #1 'repeat-authors := }
  'skip$
  if$
  current-authors 'previous-authors :=
}

ITERATE {set.current.authors}

EXECUTE {initialize.previous.authors}

ITERATE {mark.repeat.authors}

EXECUTE {initialize.previous.authors}

REVERSE {mark.repeat.authors}
%% JB added >|

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
%  "\newcommand {\unisep}{$\bullet$}" write$ newline$
  "%" write$ newline$
  "% This bibliography was produced by using jox.bst" write$ version write$ newline$
  "%" write$ newline$
  "\begin{thebibliography}{}" write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
