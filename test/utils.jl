using Compat

macro test_nothrow(ex)
    quote
        @test begin
            $(esc(ex))
            true
        end
    end
end

function nullshow(x)
    show(devnull, x)
end

function nullshow(mime::MIME, x)
    show(devnull, mime, x)
end
