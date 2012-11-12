package colormap;

public interface Mapper<I, O>
{
	public O map(I input);
}
