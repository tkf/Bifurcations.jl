macro test_nothrow(ex)
    quote
        @test begin
            $(esc(ex))
            true
        end
    end
end
