using Documenter, AstroLib

include("pages.jl")

DocMeta.setdocmeta!(AstroLib, :DocTestSetup, :(using AstroLib), recursive=true)

makedocs(;
    modules = [AstroLib],
    sitename = "AstroLib",
    format = Documenter.HTML(
        size_threshold = 400 * 1024, # 400 KiB because API reference page is big
        size_threshold_warn = 200 * 1024,
        assets = [
            "assets/favicon.ico",
        ],
    ),
    pages,
)

deploydocs(;
    repo = "github.com/JuliaAstro/AstroLib.jl.git",
    push_preview = true,
)
