# Ergänze Tabellen mit Personendaten um Angaben aus der Gemeinsamen Normdatei
# der Deutschen Nationalbibliothek
#
# Projektseite zur Rechercheoberfläche
# https://lobid.org/gnd
#
# und GUI http://lobid.org/gnd/search
#
# Beschreibung des API
# http://lobid.org/gnd/api
#
# API
# http://lobid.org/gnd/search
#
# Dokumentation
# /../README.md
#

# TODO kopiere die Dateien, git commit, git branch dev, git checkout dev
# TODO ersetze gatz mit qd
# TODO prüfe die 

"""
    module GND

Finde die GND-ID und weitere Daten zu den Bischöfen in
Gatz: Die Bischöfe des Heiligen Römischen Reiches 1198 bis 1448
"""
module GND

using DataFrames
using HTTP
using JSON
using Infiltrator
using Dates
using Logging

include("./Util.jl")


# Suche mit einer allgemeinen Suchfunktion noch Kombinationen von Feldern. Wenn ein
# Feld leer ist, wird es nie in die Suhe einbezogen.
# Sammle alle Treffer.
#
# Prüfe die einzelnen Felder. Lege ein Sortierreihenfolge fest.


"""
    QRecord

Element im Verzeichnis "member" in einer Anfrage an die Personendatenbank.
"""
struct QRecord
    data::Dict{String, Any}
end

const TOLERANCELEFTSDA = 0 # Toleranz für Amtsende(Gatz) - Sterbedatum(GND).

const TUNINGSCORE = QRecord(Dict("muster" => "fn vn sda ba", "zeitguete" => 0))
const GNDSTRINGCOLS = [:ID_GND_GND,
                       :pName_GND,
                       :Amt_GND,
                       :Sterbedatum_GND,
                       :Typ_GND,
                       :Abfrage_GND,
                       :Qualitaet_GND,
                       :nTreffer_GND,
                       :ID_VIAF_GND,
                       :ID_Wikidata_GND,
                       :URL_Wikipedia_GND]

const RANKMAX = 199

"""
    "Falsche" Treffer
"""
const SKIPIDS = [3815, #1043-1048 Eberhard von Augsburg 133613275 fn vn sda ba ip
           4635, #1137-1141 Udo (von Steinfurt?) 136876471 fn vn sda ba ip
                 ]

# Globale Variablen
"""
    COLNAMEID

Spaltenname der ID
"""
COLNAMEID = :ID_Bischof

"""
    occupations

Liste gültiger Ämter
"""
global occupations = ["Bischof", "Elekt", "Administrator", "Patriarch", "Metropolit"]
global rgxocc = Regex(join(occupations, "|"), "i")

function setoccupations(occ::Vector{<:AbstractString})
    global occupations = occ
    global rgxocc = Regex(join(occupations, "|"), "i")
end

setcolnameid(id::Symbol) = (global COLNAMEID = id)

global minscore = QRecord(Dict("muster" => "vn sda ao", "zeitguete" => 0))
#const minscore = QRecord(Dict("muster" => "", "zeitguete" => 0))


"""
    logpath

Wenn der Pfad nicht leer ist, werden hierhin die Log-Mitteilungen geschrieben.
"""
global logpath = ""
function setlogpath(logfile::AbstractString)
    path = splitdir(logfile)[1]
    path == "" && (path = ".")
    if ispath(path)
        global logpath = logfile
    else
        @error "Pfad nicht gefunden"
    end
end

global filelog = Logging.NullLogger()
global finfiltrate = false


import Base.getindex
getindex(r::QRecord, key::AbstractString) = getindex(r.data, key)

import Base.haskey
haskey(r::QRecord, key) = haskey(r.data, key)

function startlog(startid)
    global iolog
    iolog = open(LOGFILE, "a")
    println(iolog, now())
    println(iolog, "Beginne mit: ", startid)
end


"""
        makedrank()

Erstelle ein Ranking für die Qualität der Treffer. Das Verzeichnis wird für
Vergleiche genutzt.
"""
function makedrank()
    keys = String[]
    for mkn in ["fn vn sd",
                "fn vn sda",
                "fn sd",
                "fn sda",
                "vn sd",
                "vn sda",
                "fn vn",
                "fn",
                "vn"],
        mkocc in [true, false],
        mkp in [true, false]

        smk = mkn
        mkocc && (smk *= " ba")
        mkp && (smk *= " ao")

        push!(keys, smk)
    end
    push!(keys, "sda ba ao")

    ltf(a, b) = length(split(a)) < length(split(b))
    sort!(keys, alg = InsertionSort, lt = ltf, rev = true)
    keyips = fill("", 2 * length(keys))
    keyips[1:2:end] .= map(k -> k * " ip", keys)
    keyips[2:2:end] .= keys
    push!(keyips, "")

    Dict(k => i for (i, k) in enumerate(keyips))
end


import Base.isless

let drank = makedrank()
    global isless
    function isless(a::QRecord, b::QRecord)
        (islessinset(a.data["muster"],  b.data["muster"], drank)
         || (a.data["muster"] == b.data["muster"]
             && a.data["zeitguete"] < b.data["zeitguete"]))
    end
end

"""
    islessinset(a, b, drank)

Vergleiche Muster nach ihrer Position in `drank`
"""
function islessinset(a, b, drank)
    pa = get(drank, a, 0)
    pb = get(drank, b, 0)
    pa == 0 && @error ("Ungültiger Schlüssel: " * a)
    pb == 0 && @error ("Ungültiger Schlüssel: " * b)
    return pb < pa
end

let drank = makedrank()
    global getrank
    getrank(key) = drank[key]

    global isvalid
    isvalid(key) = haskey(drank, key)
end

"""
    getGND(searchtext)

Frage GND ab.

## Rückgabewert
Verzeichnis mit den Schlüsseln:
    "@context"
    "id"
    "totalItems"
    "member"
    "aggregation"
"""
function getGND(searchtext)
    st = encodeurl(searchtext)
    rq = HTTP.request("GET", "http://lobid.org/gnd/search?q=$(st)", readtimeout = 30, retries = 5);
    rq.status == 200 || error("Abfrage für '$searchtext' fehlgeschlagen")
    rds = String(rq.body);
    rds4parser = replace(rds, r"\n *" => "");
    rdt = JSON.parse(rds4parser)

end

"""
    getGND(row::DataFrameRow, col)

Frage GND ab.

## Rückgabewert
Verzeichnis mit den Schlüsseln:
    "@context"
    "id"
    "totalItems"
    "member"
    "aggregation"
"""
function getGND(row::DataFrameRow, colset)
    rdt = Dict{String, Any}()
    # Stelle keine leere Anfrage
    all(isequal(""), row[colset]) && return rdt
    searchtext = join(row[colset], " ")
    st = encodeurl(searchtext)
    rq = HTTP.request("GET", "http://lobid.org/gnd/search?q=$(st)", readtimeout = 30, retries = 5);
    rq.status == 200 || error("Abfrage für '$searchtext' fehlgeschlagen")
    rds = String(rq.body);
    rds4parser = replace(rds, r"\n *" => "");
    rdt = JSON.parse(rds4parser)
    # Schreibe die Abfrage mit in die Ergebnisse
    for rec in rdt["member"]
       rec["abfrage"] = join(string.(colset), ", ")
    end
    rdt
end


"""
    reconcile!(df::AbstractDataFrame,
                    qcols,
                    dfocc::AbstractDataFrame;
                    colswitch = :VOID,
                    nmsg = 40,
                    toldateofdeath = 5,
                    toloccdeath = 5)

Frage GND nacheinander nach den Werten aus den Spaltenkombinationen in `qcols`
ab. Verwende `toldateofdeath` als Toleranz für das Sterbedatum und
`toloccdeath` als Toleranz zwischen Amtsende und Sterbedatum. Bearbeite nur
Datensätze, für die `colswitch` nicht befüllt ist.

Beispiel für `qcols`
```
qcols = [[:Familienname, :Vorname], [:Famlienname, :Amt]]
```
"""
function reconcile!(df::AbstractDataFrame,
                    qcols,
                    dfocc::AbstractDataFrame;
                    colswitch = :VOID,
                    nmsg = 40,
                    toldateofdeath = 5,
                    toloccdeath = 5)

    global logpath
    global occupations
    startid = df[1, COLNAMEID]
    minkey = minscore["muster"]
    minrank = getrank(minkey)


    if logpath != ""
        logio = open(logpath, "a")
        filelog = SimpleLogger(logio)
        with_logger(filelog) do
            @info now() startid
            @info "Level:" minkey minrank
            @info "Gültige Ämter:" occupations
        end
    end

    @info "Level:" minkey minrank
    @info "Gültige Ämter" occupations
    

    function writematch!(row, bestrec, records)
        writerow!(row, bestrec)
        nbest = count(isequal(bestrec), records)
        if nbest > 1
            with_logger(filelog) do
                @info "Mehrere Treffer:" row[COLNAMEID] bestrec["muster"]
            end
        end
        nbest
    end

    fnotempty(v) = !ismissing(v) && !isequal(v, "")

    dbest = Dict{Int, Int}()
    irow = 1
    for row in eachrow(df)
        records = QRecord[]
        irow % nmsg == 0 && println("Datensatz: ", irow)
        irow += 1

        iep = row[COLNAMEID]
        # rette Daten
        colswitch != :VOID && row[colswitch] != "" && continue
        # überspringe Datensätze, die bekanntermaßen falsche Treffer liefern
        iep in SKIPIDS && continue

        # Ämter für row
        places = String[]
        deocs = String[]
        dflocc = Util.rowselect(dfocc, iep, COLNAMEID)
        if size(dflocc, 1) == 0
            @info "Keine Amtsdaten gefunden" iep
        else
            places = filter(fnotempty, dflocc[!, :Bistum])
            deocs = filter(fnotempty, dflocc[!, :Amtsende])
            if length(deocs) == 0
                @info "Keine Amtsendedaten gefunden" iep
            end
        end

        ffound = false
        nbest = 0
        for cols in qcols
            nbest = 0
            qres = getGND(row, cols)
            length(qres) == 0 && continue

            append!(records, evaluate!(qres, row, places, deocs, toldateofdeath, toloccdeath))

            if length(records) > 0
                bestrec, posbest = findmax(records)
                if bestrec >= TUNINGSCORE
                    nbest = writematch!(row, bestrec, records)
                    ffound = true
                    break
                end
            end
        end

        if !ffound && length(records) > 0
            bestrec, posbest = findmax(records)
            if bestrec >= minscore
                nbest = writematch!(row, bestrec, records)
            end
        end
        dbest[nbest] = get(dbest, nbest, 0) + 1
    end

    if logpath != ""
        flush(logio)
        close(logio)
    end

    dbest
end

function selectmatches(row, records)
    nbest = 0
    ffound = false
    length(records) == 0 && return nbest, ffound

    bestrec, posbest = findmax(records)
    if bestrec >= minscore
        # schreibe einen Datensatz, auch wenn er nicht eindeutig ist.
        writerow!(row, bestrec)
        nbest = count(isequal(bestrec), records)
        if nbest == 1
            ffound = true
        else
            println("Mehrere Treffer für :", row[COLNAMEID], " mit: ", bestrec["muster"])
        end
    end
    return nbest, ffound
end

"""
    evaluate!(qres::Dict{String, Any}, row, places, toldateofdeath, toloccdeath)

Bewerte die Datensätze in `qres`.
"""
function evaluate!(qres::Dict{String, Any},
                   row,
                   places,
                   deocs,
                   toldateofdeath,
                   toloccdeath)
    idep = row[COLNAMEID]
    records = QRecord[]
    nvalid = 0
    matchkey = ""
    record = Dict{String, Any}()
    sreject = String[]
    for record in qres["member"]
        # println(record["preferredName"], ": mk: ", mk, " scdod: ", scdod)

        # Die Bewertung wird in `record` in den Feldern `muster` und `zeitguete`
        # abgelegt.
        record["amuster"] = String[]
        record["zeitguete"] = 0
        isperson(record) || continue

        evaluategnfn!(record, row)
        evaluatedod!(record, idep, row[:Sterbedatum], deocs, toldateofdeath, toloccdeath)
        if record["zeitguete"] == -1
            push!(sreject, join(record["amuster"], " "))
            continue
        end
        evaluateocc!(record) # Kommt eines der Ämter vor?!
        evaluateplace!(record, places)
        evaluateip!(record)

        matchkey = record["muster"] = join(record["amuster"], " ")
        if isvalid(matchkey)
            push!(records, QRecord(record))
        else
            matchkey != "" && push!(sreject, matchkey)
        end
    end
    if length(sreject) > 0
        with_logger(filelog) do
            @info "verworfen:" idep, unique(sreject)
        end
    end
    return records
end

function isperson(record)
    reduce(|, (map(val -> occursin("Person", val), record["type"])))
end


"""
    queryGND!(row, qcols, toldateofdeath, tolendofoccupation)

Verwende `toldod` für die Bewertung des Wertes in `:Sterbedatum` und
`toleoo` für die Bewertung des Wertes in `:Amtsende`.

Es ist möglich, dem Datensatz auch Felder für die Suche hinzuzufügen
und sie immer mit dem gleichen Wert zu Belegen, z.B. "Bischof".
"""
function queryGND!(row, qcols, toldateofdeath, tolendofoccupation)
    records = Any[]
    mkeys = String[]
    searchtext = ""
    nvalid = 0
    mkey = ""
    record = nothing
    for colset in qcols
        all(isequal(""), row[colset]) && continue  # Suche nicht nach leeren Feldern
        searchtext = join(row[colset], " ")
        # println(searchtext)
        gndres = getGND(searchtext)
        record, mkey, nvalid = evaluate(gndres, row, toldateofdeath, tolendofoccupation)
        # Brich ab, sobald es gültige Treffer gibt. Dazu dürfen die Anforderungen an einen
        # gültigen Treffer nicht zu gering sein! siehe minscore
        (nvalid > 0) && break
    end
    if nvalid > 0
        writerow!(row, record, mkey, nvalid)
    end
end

"""
Bewerte die Datensätze in `gndres` anhand der Werte in `row`
"""
function evaluate(gndres, row, toldod, toleoo)
    acriteria = Pair{Matchkey, Int}[]
    nvalid = 0
    record = Dict{String, Any}()
    arec = typeof(record)[]
    amk = String[]
    mkout = ""
    mk = ""
    sdsc = 0 # Bewertung der Abweichung für das Sterbedatum
    for rcel in gndres["member"]
        amk = String[]
        mk = ""
        # println(rcel["preferredName"], ": mk: ", mk, " scdod: ", scdod)
        mk, scdod = evaluatedod(rcel, row, toldod, toleoo)
        (mk != "") || continue # Ein passendes Datum ist Grundvoraussetzung
        push!(amk, mk)
        push!(amk, evaluategnfn(rcel, row))
        push!(amk, evaluateocc(rcel))
        push!(amk, evaluateplace(rcel, row))

        mk = join(filter(!isequal(""), amk), " ")
        push!(acriteria, Pair(Matchkey(mk), scdod))
        push!(arec, rcel)
    end
    nvalid = length(acriteria)
    # ## TODO 2020-01-29 in aufrufende Funktion verlegen
    if nvalid > 0
        maxval, indmax = findmax(acriteria)
        if maxval > MINSCORE
            record = arec[indmax]
            mkout = maxval[1]
        else
            nvalid = 0
        end
    end
    return (record, mkout, nvalid)
end

const RYEAR = "([0-9]?[0-9]?[0-9]{2})(/[0-9]?[0-9]\\??)?"
const RGXYEAR = Regex(RYEAR)
const RGXDOBAD = Regex("(" * RYEAR * ")"
                     * " *((-|bis)[^1-9]*"
                     * "(" * RYEAR * "))?")


"""

Verwende `toldod` für die Bewertung der Sterbedatum in `dod` und
`toleoo` für die Bewertung der Werte in `deocs` (Amtsendedaten)

"""
function evaluatedod!(rcel,
                      id,
                      dod::Union{AbstractString, Missing},
                      deocs,
                      toldod,
                      toleoo)
    # touch 2020-02-04
    # Für die Bischöfe vor 1198 gibt es oft nur eine Angabe für das
    # Jahrhundert. "[4. Jh.]"


    function matchdod(dod,
                      datecand::Int,
                      matchkey::AbstractString,
                      toleranceleft::Int,
                      toleranceright::Int,)

        if ismissing(dod) || datecand == 0
            return 0, ""
        end

        delta = dod - datecand
        if delta <= toleranceleft && -delta <= toleranceright
            score = delta > 0 ? toleranceleft - delta : toleranceright + delta
            return score, matchkey
        else
            return -1, "dsd"
        end

    end

    # Lies Daten und bestimme maximale Qualität der Abfrage

    datecand = 0
    if haskey(rcel, "dateOfDeath")
        sdatecand = rcel["dateOfDeath"][1]
        rgm = match(RGXYEAR, sdatecand)
        if rgm != nothing
            datecand = parse(Int, rgm[1])
        else
            @warn ("Kein Datum gefunden in '" * sdatecand * "' für " * string(id))
            datecand = 0
        end
    elseif haskey(rcel, "dateOfBirthAndDeath")
        sdatecand = rcel["dateOfBirthAndDeath"][1]
        rgm = match(RGXDOBAD, sdatecand)
        if rgm != nothing
            if rgm[7] != nothing
                datecand = parse(Int, rgm[7])
            elseif rgm[2] != nothing
                datecand = parse(Int, rgm[2])
            end
        else
            datecand = 0
        end
    end

    matchkey = ""
    score = 0
    
    sdates = String[]
    if !ismissing(dod) && dod != ""
        sdates = [dod]
        toleranceleft = toleranceright = toldod
        matchkey = "sd"
    elseif length(deocs) > 0
        sdates = deocs
        toleranceleft = TOLERANCELEFTSDA
        toleranceright = toleoo;
        matchkey = "sda"
    else # keine sinnvollen Vergleichsdaten
        return matchkey, score
    end    

    # Vergleiche; finde das Maximum
    scmks = Tuple{Int, String}[]
    for sdate in sdates
        dod = parsedate(sdate)
        score, matchkey = matchdod(dod, datecand, matchkey, toleranceleft, toleranceright)
        push!(scmks, (score, matchkey))
    end

    if length(scmks) > 0
        score, matchkey = maximum(scmks)
    end

    if matchkey != ""
        push!(rcel["amuster"], matchkey)
        rcel["zeitguete"] = score
    end

    return matchkey, score
end


"""

Suche nach Familienname und oder Vorname
"""
function evaluategnfn(rcel, row)::String
    ffn = false
    fgn = false
    # Falls es keinen Familiennamen gibt, wird von `checkkey` zurückgegeben
    familyname = row[:Familienname]
    givenname = row[:Vorname]
    if haskey(rcel, "preferredNameEntityForThePerson")
        gd = rcel["preferredNameEntityForThePerson"]
        ffn |= (checkkey(gd, "nameAddition", familyname)
                || checkkey(gd, "surname", familyname))
        fgn |= (checkkey(gd, "forname", givenname)
                || checkkey(gd, "personalName", givenname))
        if ffn && fgn
            return "fn vn"
        end
    end

    if haskey(rcel, "variantNameEntityForThePerson")
        for var in rcel["variantNameEntityForThePerson"]
            ffn |= (checkkey(var, "nameAddition", familyname)
                   || checkkey(var, "surname", familyname))
            fgn |= (checkkey(gd, "forname", givenname)
                    || checkkey(var, "personalName", givenname))
            if ffn && fgn
                return "fn vn"
            end
        end
    end

    if haskey(rcel, "preferredName")
        pn = rcel["preferredName"]
        ffn |= (familyname != "" && occursin(familyname, pn))
        fgn |= occursin(givenname, pn)
        if ffn && fgn
            return "fn vn"
        end
    end

    if haskey(rcel, "variantName")
        ffn |= checkkey(rcel, "variantName", familyname)
        fgn |= checkkey(rcel, "variantName", givenname)
        if ffn && fgn
            return "fn vn"
        end
    end

    ffn && return "fn"
    fgn && return "vn"
    return ""
end

"""

Suche nach Familienname und oder Vorname
"""
function evaluategnfn!(rcel, row)::String

    function matchgnfn(rcel, row)
        # Falls es keinen Familiennamen gibt, wird von `checkkey` zurückgegeben
        ffn = false
        fgn = false

        familyname = row[:Familienname]
        givenname = row[:Vorname]
        if haskey(rcel, "preferredNameEntityForThePerson")
            gd = rcel["preferredNameEntityForThePerson"]
            ffn |= (checkkey(gd, "nameAddition", familyname)
                    || checkkey(gd, "surname", familyname))
            fgn |= (checkkey(gd, "forname", givenname)
                    || checkkey(gd, "personalName", givenname))
            if ffn && fgn
                return "fn vn"
            end
        end

        if haskey(rcel, "variantNameEntityForThePerson")
            for var in rcel["variantNameEntityForThePerson"]
                ffn |= (checkkey(var, "nameAddition", familyname)
                        || checkkey(var, "surname", familyname))
                fgn |= (checkkey(gd, "forname", givenname)
                        || checkkey(var, "personalName", givenname))

                if ffn && fgn
                    return "fn vn"
                end
            end
        end

        if haskey(rcel, "preferredName")
            pn = rcel["preferredName"]
            ffn |= (familyname != "" && occursin(familyname, pn))
            fgn |= occursin(givenname, pn)
            if ffn && fgn
                return "fn vn"
            end
        end

        if haskey(rcel, "variantName")
            ffn |= checkkey(rcel, "variantName", familyname)
            fgn |= checkkey(rcel, "variantName", givenname)
            if ffn && fgn
                return "fn vn"
            end
        end
        ffn && return "fn"
        fgn && return "vn"

        return ""
    end
    matchkey = matchgnfn(rcel, row)
    if matchkey != ""
        push!(rcel["amuster"], matchkey)
    end
    return matchkey
end


function evaluateocc(rcel)::String
    focc = false
    if haskey(rcel, "professionOrOccupation")
        ga = rcel["professionOrOccupation"]
        for gd in ga
            gp = gd["label"]
            focc = occursin(r"((b|B)ischof|(E|e)lekt)", gp)
            focc && return "ba"
        end
        return ""
    else
        return ""
    end
end

function evaluateocc!(rcel)::Bool
    # Es soll genügen, wenn überhaupt ein passendes
    # (z.B. bischöfliches) Amt ausgeübt wurde.
    global rgxocc

    focc = false
    if haskey(rcel, "professionOrOccupation")
        ga = rcel["professionOrOccupation"]
        for gd in ga
            gp = gd["label"]
            rgm = match(rgxocc, gp)
            focc = rgm != nothing
            focc && break
        end
    elseif haskey(rcel, "professionOrOccupationAsLiteral")
        for gp in rcel["professionOrOccupationAsLiteral"]
            rgm = match(rgxocc, gp)
            focc = rgm != nothing
            focc && break
        end
    end

    if focc
        push!(rcel["amuster"], "ba")
    end
    focc
end


function evaluateplace(rcel, row)::String
    fplace = false
    diocese = row[:Bistum]
    if diocese != ""  && haskey(rcel, "placeOfActivity")
        ga = rcel["placeOfActivity"]
        for gd in ga
            gp = gd["label"]
            focc = occursin(diocese, gp)
            focc && return "ao"
        end
        return ""
    else
        return ""
    end
end

"""
    Kommt einer der Orte in `places` in einer der Ortsangaben
    in `rcel` vor?
"""
function evaluateplace!(rcel, places)::Bool
    fplace = false

    if haskey(rcel, "placeOfActivity")
        ga = rcel["placeOfActivity"]
        for place in places, gd in ga
            gp = gd["label"]
            fplace = !ismissing(place) && occursin(place, gp)
            fplace && break
        end
    end
    if fplace
        push!(rcel["amuster"], "ao")
    end
    fplace
end


"""
    filterbynamedod(dgnd::Dict{String, Any},
                         givenname::AbstractString,
                         familyname::AbstractString,
                         dateofdeath::AbstractString;
                         tolerance::Int = 10)

Filtere die Einträge in `dgnd` gegen `givenname`, `familyname` und
`dateofdeath`. Das Sterbedatum darf nicht um mehr als `tolerance`
abweichen. Einheit: Jahre.

## Rückgabewert
Liste von Verzeichnissen, jeweils eins für einen Treffer
"""
function filterbynamedod(dgnd::Dict{String, Any},
                         givenname::AbstractString,
                         familyname::AbstractString,
                         dateofdeath::AbstractString;
                         tolerance::Int = 10)

    res = Array{Dict{String, Any}, 1}(undef, 0)

    # Wandle das Sterbedatum in eine Zahl um
    mdod = match(r"[1-9][0-9]+", dateofdeath)
    mdod != nothing || error("Ungültiges Sterbedatum '$(dateofdeath)'")
    idod = parse(Int, mdod.match)

    for mbel in dgnd["member"]
        mid = match(r"/([[:alnum:]-]+)$", mbel["id"])
        if mid == nothing
            @warn string("Keine ID gefunden für '", givenname, ", ", familyname, "'.")
            continue
        else
            id = mid[1]
        end

        if (findnames(mbel, givenname, familyname, id)
            && acceptdod(mbel, idod, tolerance))
            push!(res, mbel)
        end
    end
    res
end

# veraltet
function encodeurl_simple(s)
    s = replace(s, " " => "%20")
    s = replace(s, "ä" => "%C3%A4")
    s = replace(s, "ö" => "%C3%B6")
    s = replace(s, "ü" => "%C3%BC")

    s = replace(s, "á" => "%C3%A1")
    s = replace(s, "é" => "%C3%A9")
    s = replace(s, "ú" => "%C3%BA")

    s = replace(s, "Ä" => "%C3%84")
    s = replace(s, "Ö" => "%C3%96")
    s = replace(s, "Ü" => "%C3%9C")
    s = replace(s, "ß" => "%C3%9F")
end

function encodeurl(s)
    cs = String[]
    for c in s
        cu = codeunits(string(c))
        for b in cu
            push!(cs, "%")
            push!(cs, bytes2hex([b]))
        end
    end
    cs
    join(cs)
end

function acceptdod(mbel,
                   dod::Int,
                   tolerance)

    if !haskey(mbel, "dateOfDeath")
        return false
    else
        gnddod = mbel["dateOfDeath"][1]
        rg = findfirst(r"^[0-9]+", gnddod)
        rg != nothing || (@info("Keine Jahreszahl in $(gnddod)"); return false)
        ignddod = parse(Int, gnddod[rg])

        delta = abs(dod - ignddod)
        return (delta < tolerance)
    end
end


"""
"""
function findnames(mbel,
                   givenname::AbstractString,
                   familyname::AbstractString,
                   id)

    prefn = mbel["preferredName"]
    if (occursin(givenname, prefn)
        && occursin(familyname, prefn))
        return true
    elseif haskey(mbel, "variantNameEntityForThePerson")
        dvn = mbel["variantNameEntityForThePerson"]
        flaggn = false
        flagfn = false
        for vnel in dvn
        flaggn = (checkkey(vnel, "nameAddition", familyname)
                  || checkkey(vnel, "surname", familyname))
            flagfn = (checkkey(vnel, "forename", givenname)
                      || checkkey(vnel, "personalName", givenname))
            if flaggn && flagfn
                return true
            end
        end
        return false
    else
        return false
    end
end


"""
Prüfe ob `svalue` (Name) in einem Array `dn[key]` vorkommt. Mindenstens ein Teil von
`svalue` muss als Wort in einem der Elemente von `dn[key]` enthalten sein.
"""
function checkkey(dn, key, name::AbstractString)
    name == "" && return false

    # nutze aus, dass die Eingangsnamen keine Präfixe mehr enthalten
    srec = get(dn, key, Any[])
    for field in srec
        Util.checkname(name, field) && return true
    end

    # teile den Namen auf und entferne Ordnungszahlen
    # nameelems = splitname(name)
    # if haskey(dn, key) && length(nameelems) > 0
    #     srec = dn[key]
    #     for field in srec
    #         common = intersect(nameelems, splitname(field))
    #         length(common) > 0 && return true
    #     end
    # end
    return false
end

"""
    splitname(name::AbstractString)

Teile `name`auf und entferne leere Felder und Ordnungszahlen
"""
function splitname(name::AbstractString)
    rgx = r"[^\w]+"

    filtercard(a) = filter(p -> !occursin(r"^[IVX]+$", p), a)
    filterempty(a) = filter(p -> p != "", a)
    split(name, rgx) |> filterempty |> filtercard
end

function countmatches(df; nmax = 500, nmsg = 50)
    dfdst = DataFrame([String, String, String, String, String],
                      [:Familienname, :Vorname, :Sterbedatum, :GND_pn, :GND_dod],
                      0)
    histB = Dict{Int, Int}()
    nrec = 1
    for row in eachrow(df)
        nrec % nmsg == 0 && println("Datensatz: ", nrec)
        nrec += 1
        nrec < nmax || break

        sterbedatum = row[:Sterbedatum]
        vorname = row[:Vorname]
        familienname = row[:Familienname]
        if sterbedatum != ""
            rdt = GND.getGND(vorname * " " * familienname)
            rfilter = GND.filterbynamedod(rdt, vorname, familienname, sterbedatum)
            nr = size(rfilter, 1)
            if nr == 1
                gnd_pn = rfilter[1]["preferredName"]
                gnd_dod = rfilter[1]["dateOfDeath"][1]
                push!(dfdst, (familienname, vorname, sterbedatum, gnd_pn, gnd_dod))
            end
            lc = get!(histB, nr, 0)
            histB[nr] = lc + 1
        end
    end
    histB, dfdst
end


"""
    findstructdata(df::AbstractDataFrame, fsearch, skey; nmsg = 40, nmax = 200)

Finde Feldnamen oder Werte im JSON-Baum.

## Beispiel
``` julia
fsearch = GND.uniontypes!
GND.findstructdata(dfgv, fsearch, "type", nmax = 200)
```
"""
function findstructdata(df::AbstractDataFrame, fsearch, skey; nmsg = 40, nmax = 200)
    nf = Set(String[])
    nrec = 1
    for row in eachrow(df)
        nrec % nmsg == 0 && println("Datensatz: ", nrec)
        nrec += 1
        nrec < nmax || break

        vorname = row[:Vorname]
        familienname = row[:Familienname]
        rdt = GND.getGND(vorname * " " * familienname)
        fsearch(nf, rdt, skey)
    end
    nf
end


"""
    unionfields!(nf, rdt, key::String)

Füge `nf` hinzu: Feldbezeichnungen der Elemente im Feld mit dem
Bezeichener `key` aus dem Datensatz `rdt`.
"""
function unionfields!(nf, rdt, key::String)
    # Aufzurufen z.B. mit "variantNameEntityForThePerson"
    for mbel in rdt["member"]
        if haskey(mbel, key)
            avne = mbel[key]
            for vel in avne
                union!(nf, keys(vel))
            end
        end
    end
end

"""
    unionprofessionvalues!(nf, rdt, key::String)

Füge `nf` hinzu: Berufsbezeichungen im Feld mit dem Bezeichner `key`
aus dem Datensatz `rdt`.
"""
function  unionprofessionvalues!(nf, rdt, key::String)
    # Aufzurufen mit "professionOrOccupation"
    for mbel in rdt["member"]
        if haskey(mbel, key)
            avne = mbel[key]
            for dve in avne
                push!(nf, dve["label"])
            end
        end
    end
end

"""
    uniontypes!(nf, rdt, key::String)

Füge `nf` hinzu: Berufsbezeichungen im Feld mit dem Bezeichner `key`
aus dem Datensatz `rdt`.
"""
function  uniontypes!(nf, rdt, key::String)
    # Aufzurufen mit "type"
    for mbel in rdt["member"]
        if haskey(mbel, key)
            avne = mbel[key]
            for ve in avne
                push!(nf, ve)
            end
        end
    end
end

# Ansatz aufgegeben: Wir müssten Strukturen speichern und dann die Suche organisieren.
"""
    saveGND(df, keys, filename::AbstractString)
"""
function saveGND(df::AbstractDataFrame,
                 keys::Array{Symbol},
                 filename::AbstractString;
                 nmsg = 40, nmax = 200)
    open(filename, "w") do iodst
        nrec = 1
        for row in eachrow(df)
            nrec % nmsg == 0 && println("Datensatz: ", nrec)
            nrec += 1
            nrec < nmax || break

            searchtext = join(row[keys], " ")
            if searchtext == " "
                nrec -= 1
                continue
            end
            # println(searchtext)
            rdt = getGND(searchtext)
            println(iodst, rdt)
        end
    end
end

"""
"""
function reconcile!(df::AbstractDataFrame; nmsg = 40, nmax = 200, offset = 1, tolerance = 5)
    ru = min(size(df, 1), nmax)
    dfv = @view df[offset:ru, :]

    irow = 1
    for row in eachrow(dfv)
        irow % nmsg == 0 && println("Datensatz: ", irow)

        givenname = row[:Vorname]
        familyname = row[:Familienname]


        if givenname != ""    # Bischofseintrag
            if familyname != ""
                reconcilegnfn!(row, tolerance, ScN + ScBA)
            else
                reconcilegn!(row, tolerance, ScVN + ScSD)
            end
        end
        irow += 1
    end
end

function reconcilegnfn!(row, tolerance, level)
    givenname = row[:Vorname]
    dateofdeath = row[:Sterbedatum]
    familyname = row[:Familienname]
    searchtext = givenname * " " * familyname

    rdt = getGND(searchtext)
    nitems = rdt["totalItems"]
    if nitems > 0
        rdtms = rdt["member"]
        sctable = Pair{Int, String}[]
        for mbel in rdtms
            push!(sctable, scoremember(mbel,
                                       givenname,
                                       familyname,
                                       dateofdeath,
                                       tolerance))
        end

        maxval, indmax = findmax(sctable)
        # println("Score: ", maxval, " bei ", indmax)
        nmaxval = count(se -> isequal(maxval[1], se[1]), sctable)
        if maxval[1] >= level
            writerow!(row, rdtms[indmax], maxval[2], nmaxval)
        else
            # Datensätze vorhanden, aber kein überzeugender Inhalt
            row[:nTreffer_GND] = string(nitems)
        end
    else
        # Versuche es mit der Suche nach Vorname, weil die Familiennamen manchmal streuen
        # dann muss aber das Sterbedatum passen
        reconcilegn!(row, tolerance, ScVN + ScSD)
    end
end

function reconcilegn!(row, tolerance, level)
    givenname = row[:Vorname]
    searchtext = givenname * " " * "Bischof"
    dateofdeath = row[:Sterbedatum]


    rdt = getGND(searchtext)
    nitems = rdt["totalItems"]
    if nitems > 0
        rdtms = rdt["member"]
        flagn1 = size(rdtms, 1) == 1
        sctable = Pair{Int, String}[]
        for mbel in rdtms
            push!(sctable, scoremember(mbel,
                                       givenname,
                                       dateofdeath,
                                       tolerance))
        end
        maxval, indmax = findmax(sctable)
        # println("Score: ", maxval, " bei ", indmax)
        nmaxval = count(se -> isequal(maxval[1], se[1]), sctable)

        if maxval[1] >= level#  siehe Festlegung von LevelVN
            writerow!(row, rdtms[indmax], maxval[2], nmaxval)
        else
            # Datensätze vorhanden, aber kein überzeugender Inhalt
            row[:nTreffer_GND] = string(nitems)
        end
    end
end

"""

Bewerte einen Eintrag mit `givenname` und `familyname`.
"""
function scoremember(mbel, givenname, familyname, dateofdeath, tolerance)

    # Handelt es sich um eine Person?
    isperson = reduce(|, (map(val -> occursin("Person", val), mbel["type"])))
    if !isperson
         return Pair(0, "")
    end

    tot = 0

    # Name
    mk = matchkeyname(mbel, givenname, familyname)
    if mk == "fn vn"
        tot += ScFNVN
    elseif mk == "n"
        tot += ScN
    else
        # Name aus Vorname und Familienname nicht gefunden
        return Pair(0, "")
    end

    # Sterbedatum
    scdateod = scoredateofdeath(mbel, dateofdeath, tolerance)
    # println("Ergebnis bei Datum: ", scdateod)
    if scdateod > 0
        tot += scdateod
        mk = mk * " " * "sd"
    elseif scdateod <= 0 # Abweichung zu groß
        return Pair(0, mk)
    elseif scdateod == 0 # keine Daten vorhanden
        ## prüfe den Ort der Aktivität
    end


    # Amt
    if valueprofession(mbel) != ""
        tot += ScBA
        mk = mk * " " * "ba"
    end

    # Ein Eintrag für eine "individualisierte Person" wird vorgezogen
    if valuetype(mbel) == "DifferentiatedPerson"
        tot += ScIP
        mk = mk * " " * "ip"
    end



    return Pair(tot, mk)
end


"""

Suche nach `givenname`und `familienname` in `mbel`. Gib "fn vn", "n" oder "" zurück.
"""
function matchkeyname(mbel, givenname, familyname)::String
    if haskey(mbel, "preferredNameEntityForThePerson")
        gd = mbel["preferredNameEntityForThePerson"]
        ffn = (checkkey(gd, "nameAddition", familyname)
               || checkkey(gd, "surname", familyname))
        fgn = (checkkey(gd, "forname", givenname)
               || checkkey(gd, "personalName", givenname))
        if ffn & fgn
            return ("fn vn")
        end
    end

    if haskey(mbel, "variantNameEntityForThePerson")
        for var in mbel["variantNameEntityForThePerson"]
            ffn = (checkkey(var, "nameAddition", familyname)
                   || checkkey(var, "surname", familyname))
            fgn = (checkkey(gd, "forname", givenname)
                   || checkkey(var, "personalName", givenname))
            if ffn & fgn
                return ("fn vn")
            end
        end

    end

    if haskey(mbel, "preferredName")
        pn = mbel["preferredName"]
        if (occursin(givenname, pn) & occursin(familyname, pn))
            return ("n")
        else
            return ("")
        end
    end

    if haskey(mbel, "variantName")
        pn = mbel["variantName"]
        if (occursin(givenname, pn) & occursin(familyname, pn))
            return ("n")
        else
            return ("")
        end
    end

    @warn "Keine Namensfelder gefunden: " * givenname * " " * familyname
    return ("")
end

"""

Suche nach `givenname` in `mbel`. Gib "vn" oder "" zurück.
"""
function matchkeyname(mbel, givenname)::String
    if haskey(mbel, "preferredNameEntityForThePerson")
        gd = mbel["preferredNameEntityForThePerson"]
        fgn = (checkkey(gd, "forname", givenname)
               || checkkey(gd, "personalName", givenname))
        if fgn
            return ("vn")
        end
    end

    if haskey(mbel, "variantNameEntityForThePerson")
        for var in mbel["variantNameEntityForThePerson"]
            fgn = (checkkey(gd, "forname", givenname)
                   || checkkey(gd, "personalName", givenname))
            if fgn
                return ("vn")
            end
        end
    end

    return ("")
end


function scoredateofdeath(mbel, dateofdeath, tolerance)::Int
    if dateofdeath == ""
        return 0
    end

    rm = match(r"[1-9][0-9]+", dateofdeath)
    rm != nothing || (@info "Ungültiges Sterbedatum '$(dateofdeath)'"; return false)
    dateref = parse(Int, rm.match)

    if haskey(mbel, "dateOfDeath")
        rm = match(r"([1-9][0-9]{3})|([1-9][0-9]{2})", mbel["dateOfDeath"][1])
        date = parse(Int, rm.match)
        delta = abs(dateref - date)
        if delta <= tolerance
            return ScSD + tolerance - delta
        else
            return -1
        end
    else
        return 0
    end
    @warn "Dieser Punkt sollte nicht erreicht werden."
    return 0
end

"""

Bewerte einen Eintrag mit `givenname`.
"""
function scoremember(mbel, givenname, dateofdeath, tolerance)

    # Handelt es sich um eine Person?
    isperson = reduce(|, (map(val -> occursin("Person", val), mbel["type"])))
    if !isperson
         return Pair(0, "")
    end

    tot = 0

    # Der Vorname muss in einem Feld für den Vornamen vorkommen
    # Beispiel: Johann Thomas, Brixen, Bischof

    # Name
    mk = matchkeyname(mbel, givenname)
    if mk == "vn"
        tot += ScVN
    else
        # Name aus Vorname und Familienname nicht gefunden
        return Pair(0, "")
    end

    # Sterbedatum
    scdateod = scoredateofdeath(mbel, dateofdeath, tolerance)
    # println("Ergebnis bei Datum: ", scdateod)
    if scdateod > 0
        tot += scdateod
        mk = mk * " " * "sd"
    elseif scdateod < 0
        return Pair(0, mk) # Abweichung zu groß
    end

    # Amt
    if valueprofession(mbel) != ""
        tot += ScBA
        mk = mk * " " * "ba"
    end

    # Indiv. Person
    if valuetype(mbel) == "DifferentiatedPerson"
        tot += ScIP
        mk = mk * " " * "ip"
    end

    return Pair(tot, mk)
end

"""

Befülle eine Ergebniszeile `row` mit Werten aus `mbel`.
"""
function writerow!(row::DataFrameRow, record::QRecord)
    # rgm = match(r"/([[:alnum:]-]+)$", get(record, "id", ""))
    # row[:ID_GND_GND] = rm == nothing ? "" : rgm[1]
    row[:ID_GND_GND] = record["gndIdentifier"]
    row[:pName_GND] = record["preferredName"]
    row[:Amt_GND] = valueprofession(record)
    row[:Sterbedatum_GND] = valuedateofdeath(record)
    row[:Typ_GND] = valuetype(record)
    row[:Qualitaet_GND] = record["muster"]
    row[:QRang_GND] = getrank(record["muster"])
    row[:Abfrage_GND] = record["abfrage"]
    aid = valueids(record)
    row[:ID_Wikidata_GND] = aid["wikidata"]
    row[:ID_VIAF_GND] = aid["viaf"]
    row[:URL_Wikipedia_GND] = valuewikipedia(record)
end


function valueprofession(record::QRecord)::String
    if haskey(record, "professionOrOccupation")
        ga = record["professionOrOccupation"]
        for gd in ga
            gp = gd["label"]
            occursin(r"((b|B)ischof|(E|e)lekt)", gp) && return(gp)
        end
        return ""
    else
        return ""
    end
end

function valuedateofdeath(record::QRecord)::String
    if haskey(record, "dateOfDeath")
        return(record["dateOfDeath"][1])
    else
        return("")
    end
end

function valuetype(record)::String
    if haskey(record, "type")
        ga = record["type"]
        for ttest in ["DifferentiatedPerson", "Person", "UndifferentiatedPerson"]
            ttest in ga && return(ttest)
        end
        return ""
    else
        return ""
    end
end

function evaluateip!(rcel)::Bool
    fip = "DifferentiatedPerson" == valuetype(rcel)
    if fip
        push!(rcel["amuster"], "ip")
    end
    fip
end


function valuewikipedia(record)
    if haskey(record, "wikipedia")
        ga = record["wikipedia"]
        if haskey(ga[1], "id")
            return ga[1]["id"]
        else
            @warn "Im Element Wikipedia 'id' nicht gefunden " * record["preferredName"]
        end
    else
        return ""
    end
end


function valueids(record)
    keys = ["wikidata", "viaf"]
    did = Dict([k => "" for k in keys])
    rgx = r"[^/]+$"

    function filldid(sameAs)
        valid = sameAs["id"]
        for key in keys
            if occursin(key * ".org", valid)
                rm = match(rgx, valid)
                if rm != nothing
                    did[key] = rm.match
                else
                    @warn ("Schlüssel " * key *
                           " im Feld 'id' nicht gefunden für " * record["preferredName"])
                end
            end
        end
    end

    haskeyid(sameAs) = haskey(sameAs, "id")

    if haskey(record, "sameAs")
        foreach(filldid, filter(haskeyid, record["sameAs"]))
    end

    return did
end


function valueidsT(record)
    # Implementierung mit Elementen von Transducres
    # nicht effizient
    keys = ["wikidata", "viaf"]
    did = Dict([k => "" for k in keys])

    function filldid(p) # p[1] key, p[2] Vergleichswert
        rgm = match(r"[^/]+$", p[2])
        if rgm != nothing
            did[p[1]] = rgm.match
        else
            @warn ("Schlüssel " * key *
                   " im Feld 'id' nicht gefunden für " * record["preferredName"])
        end
    end

    function idtoval(gd)
        foreach(filldid,
                Map(key -> (key, gd["id"]))
                |> Filter(matchkey),
                keys)
    end


    matchkey(p) = occursin(p[1] * ".org", p[2])

    if haskey(record, "sameAs")
        foreach(idtoval,
                Filter(gd -> haskey(gd, "id")),
                record["sameAs"])
    end
    return did
end


function makeGNDDataFrame(df)
    # Übernehme alle Daten
    dfgnd = copy(df)

    for col in GNDSTRINGCOLS
        dfgnd[!, col] .= ""
    end
    dfgnd[!, :QRang_GND] .= RANKMAX

    dfgnd
end

"""
Zähle nicht-leere Datensätze in `df`. Wenn nur einzelne Spalten gewünscht sind,
übergebe man eine entsprechende View auf `df`.
"""
function countnotempty(df::AbstractDataFrame)
    ncol = length(names(df))
    dfc = DataFrame(fill(Int, ncol), names(df), 0)
    counts = fill(0, ncol)
    for row in eachrow(df), (i, v) in enumerate(row)
        v != "" && (counts[i] += 1)
    end
    push!(dfc, counts)
    dfc
end

"""
    clearbylevel!(df::AbstractDataFrame, matchkey::AbstractString, level::Int)

Lösche Datensätze, deren Qualität unter `matchkey`, bzw. `rank` liegt. Die Angaben
sind zur Sicherheit redundant.
"""
function clearbylevel!(df::AbstractDataFrame, matchkey::AbstractString, rank::Int)
    if rank != getrank(matchkey)
        error("Muster und Rang stimmen nicht überein")
    end
    colrank = df[!, :QRang_GND]
    ixr = colrank .> rank
    df[ixr, GNDSTRINGCOLS] .= "" # zulässig!
    df[ixr, :QRang_GND] .= RANKMAX
    count(df[!, :QRang_GND] .< RANKMAX) # Zahl der Datensätze, die erhalten blieben
end

"""
    clearbyindex!(df::AbstractDataFrame, ixr)

Lösche alle GND-Daten aus `df` mit Index `ixr`.
"""
function clearbyindex!(df::AbstractDataFrame, ixr)
    if size(df, 1) != length(ixr)
        @error Unterschiedliche Länge der Argumente
    end
    df[ixr, GNDSTRINGCOLS] .= "" # zulässig!
    df[ixr, :QRang_GND] .= RANKMAX
    count(ixr)
end


"""
    infolevel!(df::AbstractDataFrame, matchkey::AbstractString)

"""
function infolevel(df::AbstractDataFrame, matchkey::AbstractString)
    rank = getrank(matchkey)
    println("Das Muster '", matchkey, "' hat die Position: ", rank, ".")
    println("Datensätze insgesamt: ", size(df, 1))
    colrank = df[!, :QRang_GND]
    println("Rang größer als '", matchkey, "': ", count(colrank .> rank))
    rankmax = maximum(colrank)
    nvalid = count(colrank .< rankmax)
    ndel = count((colrank .< rankmax) .& (colrank .> rank))
    println("Davon kleiner maximaler Rang : ", ndel)
    println("Es bleiben erhalten: ", nvalid - ndel)
end

"""
    deltaID(ref, gs)

Vergleiche IDs: `ref` ist z.B. manuell bestätigt.

## Beispiel
```
deltaID(r[[:GND_ID, :ID_GND_GND]]...)
```
"""
function deltarow(vref, vq)
    if vq == ""
        vref == "" ? (return "vv") : (return "xv")
    elseif vref == ""
        return "vx"
    else
        vref == vq ? (return "GLEICH") : (return "UNGLEICH")
    end
end

"""
    getsda(dfgnd::AbstractDataFrame,
           dfocc::AbstractDataFrame)

Berechne die Differenz zwischen Sterbedatum und letztem Amt
"""
function getsda(dfgnd::AbstractDataFrame,
                dfocc::AbstractDataFrame)
    typeid = eltype(dfgnd[!, COLNAMEID])
    df = DataFrame([typeid, Int], [COLNAMEID, :sda])
    rowsda = Dict{Symbol, Union{typeid, Int}}()

    for row in eachrow(dfgnd)
        sdateod = row[:Sterbedatum_GND]
        iep = row[COLNAMEID]

        rowsda[COLNAMEID] = iep

        dateod = parseyear(sdateod)
        dateod == nothing && continue

        csoe = @from i in dfocc begin
            @where i.ID_Bischof == iep
            @select i.Amtsende
            @collect
        end

        length(csoe) == 0 && continue

        # Zähle nur 0 und größer
        csda = Int[]
        for soe in csoe
            dateoe = parseyear(soe)
            dateoe == nothing && continue
            #println("  Amtsende: ", dateoe)
            delta = dateod - dateoe
            delta >= 0 && push!(csda, delta)
        end
        if length(csda) > 0
            rowsda[:sda] = minimum(csda)
            push!(df, rowsda)
        end
    end
    return df
end


function parseyear(s)
    rgm = match(RGXYEAR, s)
    if rgm == nothing
        return nothing
    else
        return parse(Int, rgm[1])
    end
end

"""

Gib `missing` zurück, wenn kein gültiges Datum gefunden werden kann
"""
function parsedate(sdate::Union{AbstractString, Missing})
    rgxyear = r"[0-9]?[0-9]?[0-9]{2}"
    valdate = missing

    ismissing(sdate) && return valdate

    sdate in ("", "(?)", "?") && return valdate
    
    rgm = match(rgxyear, sdate)
    if rgm == nothing
        @warn ("Ungültiges Datum in: " * sdate)
    else
        valdate = parse(Int, rgm.match)
    end

    return valdate
end

function parsedate(sdate::Int)
    sdate
end




end # modul GND
