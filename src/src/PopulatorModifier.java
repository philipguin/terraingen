package src;

public class PopulatorModifier<T> implements IPopulator<T>
{
	private final IPopulator<? super T> populator;
	private final IModifier<? super T> modifier;
	
	public PopulatorModifier(IPopulator<? super T> populator, IModifier<? super T> modifier)
	{
		this.populator = populator;
		this.modifier = modifier;
	}
	
	public final void populate(T values)
	{
		populator.populate(values);
		modifier.modify(values);
	}
}
